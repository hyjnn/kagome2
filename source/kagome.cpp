#include <bitset>
#include <cctype> //std::isdigit
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional> //std::function
#include <type_traits>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>

#include <Eigen/Dense>

#include "Polynomial.h"
#include "State.h"
#include "Timer.h"

using FloatType = double;
using FloatVector = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
using FloatMatrix = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename _Scalar>
concept EigenScalar = requires (Eigen::Matrix<_Scalar, -1, -1> testMat1, Eigen::Matrix<_Scalar, -1, -1> testMat2) {
    testMat1 * testMat2;
};

unsigned long binomialCoefficient(unsigned n, unsigned k) {
    if (k == 0 || k == n) {
        return 1;
    }
    if (k > n) {
        return 0;
    }
    return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
}

std::tuple<FloatType, FloatVector> maxEigen(FloatMatrix mat, unsigned n, FloatVector testVect = {}) {
    if (testVect.size() == 0) testVect = FloatVector::Random(mat.cols());
    testVect /= testVect.norm();
    FloatType eigenvalue = 0;

    for (unsigned i = 0; i < n; i++) {
        testVect = mat * testVect;
        eigenvalue = testVect.norm();

        testVect /= eigenvalue;
    }

    return std::make_tuple(eigenvalue, testVect);
}

//Returns the amount of x and y type nodes respectively as fractions
std::pair<Polynomial::CoeffScalar, Polynomial::CoeffScalar> countTypesKagome(Polynomial p) {
    std::pair<Polynomial::CoeffScalar, Polynomial::CoeffScalar> result = { 0, 0 };
    Polynomial::CoeffScalar sum = 0;

    for (size_t i = 0; i <= p.degree(); i++) {
        for (size_t j = 0; j <= p.degree() - i; j++) {
            if (p.coefficients(i, j) != 0) {
                result.first += p.coefficients(i, j) * i;
                result.second += p.coefficients(i, j) * j;
                sum += p.coefficients(i, j);
            }
        }
    }
    result.first /= sum * p.degree();
    result.second /= sum * p.degree();

    return result;
}

FloatType countSqStates(const State& inState, const State& outState) {
    FloatType exchangeType = -1;
    FloatType firstExchange = -1;

    for (std::size_t k = 0; k < inState.size(); k++) {
        if (inState.at(k) != outState.at(k)) {
            if (exchangeType == -1) firstExchange = inState.at(k);
            if (exchangeType == inState.at(k)) return 0;
            exchangeType = inState.at(k);
        }
    }

    if (firstExchange == -1) return 2;
    if (firstExchange == exchangeType) return 0;
    return 1;
}

FloatType countDcStates(const State& inState, const State& outState) {
    if (inState.size() % 2 != 0)
        throw std::invalid_argument("Wymiar musi byc parzysty, aby moc obliczyc macierz dla nieodzialujacych atomow");

    if (inState.at(0) != outState.at(outState.size() - 1) || inState.at(1) != outState.at(0))
        if (inState.at(0) != outState.at(0) || inState.at(1) != outState.at(outState.size() - 1))
            return 0;
    for (std::size_t k = 2; k < inState.size(); k += 2)
        if (inState.at(k) != outState.at(k - 1) || inState.at(k + 1) != outState.at(k))
            if (inState.at(k) != outState.at(k) || inState.at(k + 1) != outState.at(k - 1))
                return 0;

    return 1;
}

//slower than square-decoupled product approach
//top receives top arrow, gets overwritten with bottom arrow
static FloatType countKgStatesInternal(const State& inState, const State& outState, FloatType& top) {
    FloatType countMultiplier = 1;

    for (size_t i = 2; i < inState.size(); i += 2) {
        if (inState.at(i) == top) {
            if (outState.at(i - 1) != outState.at(i)) {
                countMultiplier *= 2;
                top = 1 - inState.at(i + 1);
            }
            else {
                if (inState.at(i + 1) != outState.at(i)) {
                    countMultiplier *= 0;
                    break;
                }
                top = inState.at(i + 1);
            }
        }
        else {
            if (outState.at(i - 1) != outState.at(i)) {
                if (inState.at(i + 1) != top) {
                    countMultiplier *= 0;
                    break;
                }
                top = inState.at(i + 1);
            }
            else {
                if (inState.at(i) != outState.at(i)) {
                    countMultiplier *= 0;
                    break;
                }
                top = 1 - inState.at(i + 1);
            }
        }
    }

    return countMultiplier;
}

//slower than square-decoupled product approach
FloatType countKgStates(const State& inState, const State& outState) {
    std::vector<FloatType> variantCounts(2);
    std::vector<FloatType> variantTops(2, 0), variantBots(2, 0); //arrow up is 1, arrow down is 0

    if (outState.at(outState.size() - 1) != outState.at(0)) {
        if (inState.at(0) == inState.at(1)) {
            variantCounts = { 2, 0 };
            variantTops.at(0) = inState.at(0);
            variantBots.at(0) = 1 - inState.at(1);
        }
        else {
            variantCounts = { 2, 1 };
            variantTops = { inState.at(0) , 1 - inState.at(0) };
            variantBots = variantTops;
        }
    }
    else {
        if (inState.at(0) == outState.at(0)) {
            if (inState.at(1) == outState.at(0)) {
                variantCounts = { 1, 1 };
                variantTops = { 1, 0 };
                variantBots = { 1, 0 };
            }
            else {
                variantCounts = { 1, 0 };
                variantTops.at(0) = inState.at(1);
                variantBots.at(0) = 1 - inState.at(1);
            }
        }
        else {
            if (inState.at(1) == outState.at(0)) {
                variantCounts = { 1, 0 };
                variantTops.at(0) = inState.at(0);
                variantBots.at(0) = 1 - inState.at(0);
            }
            else {
                if (inState.at(0) != inState.at(1)) {
                    variantCounts = { 1, 0 };
                    variantTops.at(0) = inState.at(0);
                    variantBots.at(0) = 1 - inState.at(1);
                }
                else {
                    variantCounts = { 0, 0 };
                }
            }
        }
    }

    for (size_t i = 0; i < variantCounts.size(); i++) {
        if (variantCounts.at(i) != 0) {
            variantCounts.at(i) *= countKgStatesInternal(inState, outState, variantBots.at(i));
            if (variantTops.at(i) != variantBots.at(i)) variantCounts.at(i) = 0;
        }
    }

    return variantCounts.at(0) + variantCounts.at(1);
}

//returns number of x type nodes, y type nodes and number of configurations respectively
static std::tuple<Polynomial::CoeffScalar, Polynomial::CoeffScalar, Polynomial::CoeffScalar> countKgNodeTypesInternal(const State& inState, const State& outState, FloatType& top) {
    Polynomial::CoeffScalar xCount = 0, yCount = 0;
    Polynomial::CoeffScalar countMultiplier = 1;

    for (size_t i = 2; i < inState.size(); i += 2) {
        if (inState.at(i) == top) {
            if (outState.at(i - 1) != outState.at(i)) {
                xCount *= 2;
                yCount *= 2;
                yCount += countMultiplier * 3;

                countMultiplier *= 2;
                top = 1 - inState.at(i + 1);
            }
            else {
                if (inState.at(i + 1) != outState.at(i)) {
                    countMultiplier *= 0;
                    break;
                }
                top = inState.at(i + 1);
                xCount += countMultiplier * 2;
                if (inState.at(i) == outState.at(i - 1)) {
                    yCount += countMultiplier;
                }
            }
        }
        else {
            if (outState.at(i - 1) != outState.at(i)) {
                if (inState.at(i + 1) != top) {
                    countMultiplier *= 0;
                    break;
                }
                top = inState.at(i + 1);
                xCount += countMultiplier * 2;
                if (inState.at(i) != outState.at(i - 1)) {
                    yCount += countMultiplier;
                }
            }
            else {
                if (inState.at(i) != outState.at(i)) {
                    countMultiplier *= 0;
                    break;
                }
                top = 1 - inState.at(i + 1);
                xCount += countMultiplier * 2;
                if (inState.at(i + 1) == outState.at(i - 1)) {
                    yCount += countMultiplier;
                }
            }
        }
    }

    return std::make_tuple(xCount, yCount, countMultiplier);
}

std::tuple<Polynomial::CoeffScalar, Polynomial::CoeffScalar> countKgNodeTypes(size_t inStateIndex, size_t outStateIndex, const State& inState, const State& outState, const FloatMatrix& tMat) {
    std::vector<std::vector<Polynomial::CoeffScalar>> variantCounts;
    std::vector<FloatType> variantTops(2, 0), variantBots(2, 0); //arrow up is 1, arrow down is 0

    if (outState.at(outState.size() - 1) != outState.at(0)) {
        if (inState.at(0) == inState.at(1)) {
            variantCounts = { { 0, 1.5, 2 }, { 0, 0, 0 } }; //Two ways to fill this in, one gives y^2z, the other yz^2, together y^3z^3. I add them up immediately, this has to be accounted for later on by dividing node counts by 2.
            variantTops.at(0) = inState.at(0);
            variantBots.at(0) = 1 - inState.at(1);
        }
        else {
            variantCounts = { { 0, 1.5, 2 }, { 2, static_cast<Polynomial::CoeffScalar>(inState.at(0) == outState.at(0)), 1}}; //First variant is the same situation as the one above
            variantTops = { inState.at(0) , 1 - inState.at(0) };
            variantBots = variantTops;
        }
    }
    else {
        if (inState.at(0) == outState.at(0)) {
            if (inState.at(1) == outState.at(0)) {
                variantCounts = { { 2, 1, 1 }, { 2, 1, 1 } };
                variantTops = { 1, 0 };
                variantBots = { 1, 0 };
            }
            else {
                variantCounts = { { 2, 0, 1 }, {0, 0, 0} };
                variantTops.at(0) = inState.at(1);
                variantBots.at(0) = 1 - inState.at(1);
            }
        }
        else {
            if (inState.at(1) == outState.at(0)) {
                variantCounts = { { 2, 0, 1 }, { 0, 0, 0 } };
                variantTops.at(0) = inState.at(0);
                variantBots.at(0) = 1 - inState.at(0);
            }
            else {
                variantCounts = { { 0, 0, 0 }, { 0, 0, 0 } };
            }
        }
    }

    if (variantCounts.at(0) != variantCounts.at(1) && variantCounts.at(1).at(2) != 0) {
        for (size_t i = 0; i < variantCounts.size(); i++) {
            if (variantCounts.at(i).at(2) != 0) {
                variantCounts.at(i).at(2) *= countKgStatesInternal(inState, outState, variantBots.at(i));
                if (variantTops.at(i) != variantBots.at(i)) variantCounts.at(i).at(2) = 0;
            }
        }
    }
    else {
        for (size_t i = 0; i < variantCounts.size() && variantCounts.at(i).at(2) != 0; i++) {
            variantCounts.at(i).at(2) = tMat(inStateIndex, outStateIndex) / variantCounts.at(i).at(2);
        }
    }

    return std::make_tuple(variantCounts.at(0).at(0) * variantCounts.at(0).at(2) + variantCounts.at(1).at(0) * variantCounts.at(1).at(2),
                            variantCounts.at(0).at(1) * variantCounts.at(0).at(2) + variantCounts.at(1).at(1) * variantCounts.at(1).at(2));
}

std::tuple<Polynomial::CoeffScalar, Polynomial::CoeffScalar> eigenGetRatios(const Eigen::MatrixXPoly& tMat, const FloatVector& eigenvector, unsigned dim) {
    std::size_t min = (std::pow(2, dim) - binomialCoefficient(dim, dim / 2)) / 2, max = (std::pow(2, dim) + binomialCoefficient(dim, dim / 2)) / 2;

    Polynomial sum(tMat(0, 0).degree());

    for (size_t i = min; i < max; i++) {
        for (size_t j = min; j < max; j++) {
            sum += eigenvector(i) * tMat(i, j) * eigenvector(j);
        }
    }

    return countTypesKagome(sum);
}

std::tuple<Polynomial::CoeffScalar, Polynomial::CoeffScalar> transferGetRatios(const FloatMatrix& tMat, unsigned dim) {
    std::size_t i = (std::pow(2, dim) - binomialCoefficient(dim, dim / 2)) / 2, maxj = i, j;
    State inState(dim), outState(dim);
    std::tuple<Polynomial::CoeffScalar, Polynomial::CoeffScalar> temp;
    Polynomial::CoeffScalar xCount = 0, yCount = 0;
    Polynomial::CoeffScalar nodeCount = 0;

    inState.begin(dim / 2);

    while (true) {
        j = maxj;
        outState.begin(dim / 2);

        do {
            temp = countKgNodeTypes(i, j, inState, outState, tMat);
            xCount += std::get<0>(temp);
            yCount += std::get<1>(temp);

            j += 1;
        } while (outState.next());

        i += 1;
        if (!inState.next()) {
            maxj = j;
            break;
        }
    }
    nodeCount *= 3;

    return std::make_tuple(xCount, yCount);
}

//row number corresponds to inState, column number to outState
template<EigenScalar _Scalar>
Eigen::Matrix<_Scalar, -1, -1> genTransfer(std::function<_Scalar(const State&, const State&)> countFunc, unsigned dim) {
    Eigen::Matrix<_Scalar, -1, -1> resultMatrix = Eigen::Matrix<_Scalar, -1, -1>::Zero(static_cast<int>(std::pow(2, dim)), static_cast<int>(std::pow(2, dim)));
    std::size_t i = 0, maxj = 0, j;
    State inState(dim), outState(dim);

    for (std::size_t m = 0; m < dim + 1; m++) {
        inState.begin(m);

        while (true) {
            j = maxj;
            outState.begin(m);

            do {
                resultMatrix(i, j) = countFunc(inState, outState);
                j += 1;
            } while (outState.next());

            i += 1;
            if (!inState.next()) {
                maxj = j;
                break;
            }
        }
    }

    return resultMatrix;
}

Polynomial countSqPoly(const State& inState, const State& outState) {
    int exchangeIndex = -1;
    bool statesIn;
    size_t xCount = 0, yCount = 0;
    
    for (size_t i = 0; i < inState.size() && exchangeIndex == -1; i++) {
        if (inState.at(i) != outState.at(i)) {
            exchangeIndex = i;
            statesIn = inState.at(i);
        }
    }

    if (exchangeIndex == -1) {
        Polynomial result(inState.size());
        result.coefficients(inState.size(), 0) = 1;
        result.coefficients(0, inState.size()) = 1;

        return result;
    }

    for (size_t i = 0; i < exchangeIndex; i++) {
        (inState.at(i) == statesIn ? xCount : yCount)++;
    }

    for (size_t i = exchangeIndex + 1; i < inState.size(); i++) {
        if (inState.at(i) == outState.at(i)) {
            (inState.at(i) == inState.at(exchangeIndex) ? yCount : xCount)++;
        }
        else {
            if (inState.at(i) == inState.at(exchangeIndex)) return Polynomial(0);

            exchangeIndex = i;
        }
    }

    if (inState.at(exchangeIndex) == statesIn) return Polynomial(0);

    Polynomial result(inState.size());
    result.coefficients(xCount, yCount) = 1;
    
    return result;
}

Polynomial countDcPoly(const State& inState, const State& outState) {
    if (inState.size() % 2 != 0)
        throw std::invalid_argument("Wymiar musi byc parzysty, aby moc obliczyc macierz dla nieodzialujacych atomow");

    size_t xCount = 0, yCount = 0;

    if (inState.at(0) == outState.at(outState.size() - 1)) {
        if (inState.at(1) == outState.at(0)) {
            if (inState.at(0) == inState.at(1)) {
                xCount++;
            }
            //if we get here its a z-type node
        }
        else {
            return Polynomial(0);
        }
    }
    else if (inState.at(0) == outState.at(0) && inState.at(1) == outState.at(outState.size() - 1)) {
        yCount++;
    }
    else {
        return Polynomial(0);
    }

    for (size_t i = 2; i < inState.size(); i += 2) {
        if (inState.at(i) == outState.at(i - 1)) {
            if (inState.at(i + 1) == outState.at(i)) {
                if (inState.at(i) == inState.at(i + 1)) {
                    xCount++;
                }
                //if we get here its a z-type node
            }
            else {
                return Polynomial(0);
            }
        }
        else if (inState.at(i) == outState.at(i) && inState.at(i + 1) == outState.at(i - 1)) {
            yCount++;
        }
        else {
            return Polynomial(0);
        }
    }

    Polynomial result(inState.size() / 2);
    result.coefficients(xCount, yCount) = 1;

    return result;
}

//helper function to use in power
template<EigenScalar _Scalar, int _Rows, int _Cols>
static Eigen::Matrix<_Scalar, _Rows, _Cols> specialMult(const Eigen::Matrix<_Scalar, _Rows, _Cols>& A, const Eigen::Matrix<_Scalar, _Rows, _Cols>& B) {
    if constexpr (std::is_same<_Scalar, Polynomial>::value) {
        return A.lazyProduct(B);
    }
    else {
        return A * B;
    }
}

template<EigenScalar _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols>& power(Eigen::Matrix<_Scalar, _Rows, _Cols>& mat, unsigned n) {
    Eigen::Matrix<_Scalar, _Rows, _Cols> original = mat, temp;

    for (unsigned i = static_cast<unsigned>(std::log2(n)); i > 0; i--) {
        mat = specialMult(mat, mat);
        if (n & (1 << i - 1)) {
            mat = specialMult(mat, original);
        }
    }

    return mat;
}

void runCalcTransfer(std::string type, unsigned nrows, unsigned ncols) {
    Timer stopwatch;
    unsigned nodeCount;
    FloatType stateCount;

    if (type.substr(type.size() - 4) != "Poly") {
        FloatMatrix tMat;

        if (type == "square") {
            stopwatch.start();
            tMat = genTransfer(std::function(countSqStates), nrows);
            stopwatch.end();
            nodeCount = nrows * ncols;
        }
        else if (type == "decoupled") {
            stopwatch.start();
            tMat = genTransfer(std::function(countDcStates), nrows);
            stopwatch.end();
            nodeCount = nrows / 2 * ncols;
        }
        else if (type == "kagome") {
            stopwatch.start();
            FloatMatrix square = genTransfer(std::function(countSqStates), nrows);
            stopwatch.end();
            std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            FloatMatrix decoupled = genTransfer(std::function(countDcStates), nrows);
            stopwatch.end();
            std::cout << std::format("Decoupled: czas liczenia : {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            tMat = square * decoupled;
            stopwatch.end();
            nodeCount = nrows / 2 * 3 * ncols;
        }
        else if (type == "kagomeAlt") {
            stopwatch.start();
            tMat = genTransfer(std::function(countKgStates), nrows);
            stopwatch.end();
            nodeCount = nrows / 2 * 3 * ncols;
        }
        else {
            throw std::invalid_argument("Invalid grid type");
        }
        std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

        if (ncols > 1) {
            stopwatch.start();
            power(tMat, ncols);
            stopwatch.end();
            std::cout << std::format("Czas potegowania: {} ms\n", stopwatch.read<Timer::ms>());
        }

        stateCount = tMat.trace();
        std::cout << "Liczba stanow: " << stateCount << '\n';
        std::cout << std::format("Entropia resztkowa na czasteczke: {} k_B\n", std::log(stateCount) / nodeCount);
        
        if (nrows < 7 && ncols < 7) {
            std::cout << tMat;
        }
    }
    else {
        Eigen::MatrixXPoly polyMat;

        if (type == "kagomePoly") {
            stopwatch.start();
            Eigen::MatrixXPoly square = genTransfer(std::function(countSqPoly), nrows);
            stopwatch.end();
            std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            Eigen::MatrixXPoly decoupled = genTransfer(std::function(countDcPoly), nrows);
            stopwatch.end();
            std::cout << std::format("Decoupled: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            polyMat = square.lazyProduct(decoupled); //regular product doesn't work for big matrices
            stopwatch.end();
        }
        else {
            throw std::invalid_argument("Invalid grid type");
        }
        std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

        if (ncols > 1) {
            stopwatch.start();
            power(polyMat, ncols);
            stopwatch.end();
            std::cout << std::format("Czas potegowania: {} ms\n", stopwatch.read<Timer::ms>());
        }

        Polynomial trace = polyMat.trace();
        if (ncols * nrows < 10 || ncols == 1) std::cout << "Wielomian: " << trace << '\n';
        auto [xFraction, yFraction] = countTypesKagome(trace);

        std::cout << std::format("Ulamek wierzcholkow typu x: {}\nUlamek wierzcholkow typu y : {}", xFraction, yFraction);
    }

    std::cout << '\n';
}

void runCalcEigen(std::string type, unsigned nrows, unsigned powerIterCount = 10) {
    Timer stopwatch;
    FloatType eigenvalue;
    unsigned nodeFactor = nrows;
    FloatVector eigenvector;
    FloatMatrix tMat;

    if (type == "square") {
        stopwatch.start();
        tMat = genTransfer(std::function(countSqStates), nrows);
        stopwatch.end();
    }
    else if (type == "decoupled") {
        stopwatch.start();
        tMat = genTransfer(std::function(countDcStates), nrows);
        stopwatch.end();
    }
    else if (type == "kagome" || type == "kagomePoly") {
        stopwatch.start();
        FloatMatrix square = genTransfer(std::function(countSqStates), nrows);
        stopwatch.end();
        std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        FloatMatrix decoupled = genTransfer(std::function(countDcStates), nrows);
        stopwatch.end();
        std::cout << std::format("Decoupled: czas liczenia : {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        tMat = square * decoupled;
        stopwatch.end();
        nodeFactor = nodeFactor * 3 / 2;
    }
    else {
        throw std::invalid_argument("Invalid grid type");
    }
    std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

    stopwatch.start();
    std::tie(eigenvalue, eigenvector) = maxEigen(tMat, powerIterCount);
    stopwatch.end();
    std::cout << std::format("Czas liczenia wartosci wlasnej: {} ms\n", stopwatch.read<Timer::ms>());
    eigenvector /= eigenvector.sum();

    std::cout << "Maksymalna wartosc wlasna: " << eigenvalue << '\n';
    std::cout << std::format("Graniczna entropia resztkowa na czasteczke: {} k_B\n", std::log(eigenvalue) / nodeFactor);

    if (nrows < 7) std::cout << "Wektor wlasny:\n" << eigenvector << '\n';

    if (type == "kagomePoly") {
        stopwatch.start();
        Eigen::MatrixXPoly square = genTransfer(std::function(countSqPoly), nrows);
        stopwatch.end();
        std::cout << std::format("Square (Poly): czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        Eigen::MatrixXPoly decoupled = genTransfer(std::function(countDcPoly), nrows);
        stopwatch.end();
        std::cout << std::format("Decoupled (Poly): czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        Eigen::MatrixXPoly polyMat = square.lazyProduct(decoupled); //regular product doesn't work for big matrices
        stopwatch.end();
        std::cout << std::format("Czas liczenia macierzy transferu (Poly): {} ms\n", stopwatch.read<Timer::ms>());

        auto [xCount, yCount] = eigenGetRatios(polyMat, eigenvector, nrows);
        std::cout << std::format("Ulamek wierzcholkow x: {}\nUlamek wierzcholkow y: {}\n", xCount, yCount);
    }
}

void runPolyToFile(unsigned nrows, unsigned startColumns, unsigned maxColumns, unsigned columnStep) {
    Timer stopwatch;
    Polynomial trace;

    stopwatch.start();
    Eigen::MatrixXPoly square = genTransfer(std::function(countSqPoly), nrows);
    stopwatch.end();
    std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

    stopwatch.start();
    Eigen::MatrixXPoly decoupled = genTransfer(std::function(countDcPoly), nrows);
    stopwatch.end();
    std::cout << std::format("Decoupled: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

    stopwatch.start();
    Eigen::MatrixXPoly polyMat = square.lazyProduct(decoupled); //regular product doesn't work for big matrices
    stopwatch.end();
    std::cout << std::format("Czas liczenia macierzy transferu {} x 1: {} ms\n", nrows, stopwatch.read<Timer::ms>());

    Eigen::MatrixXPoly stepMat = polyMat;

    if (startColumns > 1) {
        stopwatch.start();
        power(polyMat, startColumns);
        stopwatch.end();
        std::cout << std::format("Czas potegowania do macierzy {} x {}: {} ms\n", nrows, startColumns, stopwatch.read<Timer::ms>());
    }

    if (columnStep > 1) {
        stopwatch.start();
        power(stepMat, columnStep);
        stopwatch.end();
        std::cout << std::format("Czas potegowania macierzy mnozacej: {} ms\n", stopwatch.read<Timer::ms>());
    }

    if (!std::filesystem::exists("./output")) std::filesystem::create_directory("./output");
    std::ofstream out(std::format("output/kagomepoly_{}x{}-{}_s{}.txt", nrows, startColumns, maxColumns, columnStep));
    trace = polyMat.trace();
    out << std::format("Siec {} x {}:\n", nrows, startColumns) << trace;

    stopwatch.start();
    auto [xFraction, yFraction] = countTypesKagome(trace);
    stopwatch.end();
    std::cout << std::format("Czas liczenia stosunkow: {}\n", stopwatch.read<Timer::ms>());
    out << std::format("\n{} {}", xFraction, yFraction);

    for (unsigned i = 1; startColumns + i * columnStep <= maxColumns; i++) {
        stopwatch.start();
        polyMat = specialMult(polyMat, stepMat);
        stopwatch.end();
        std::cout << std::format("Czas liczenia macierzy {} x {}: {} ms\n", nrows, startColumns + i * columnStep, stopwatch.read<Timer::ms>());
        trace = polyMat.trace();
        out << std::format("\nSiec {} x {}:\n", nrows, startColumns + i * columnStep) << trace;

        stopwatch.start();
        auto [xFraction, yFraction] = countTypesKagome(trace);
        stopwatch.end();
        std::cout << std::format("Czas liczenia stosunkow: {}\n", stopwatch.read<Timer::ms>());
        out << std::format("\n{} {}", xFraction, yFraction);
    }
}

std::vector<std::string> split(std::string str) {
    std::string word;
    std::vector<std::string> result;
    size_t start = 0, end = str.find(' ');

    while (end != std::string::npos) {
        word = str.substr(start, end - start);

        if (word != "") {
            result.push_back(word);
        }

        start = end + 1;
        end = str.find(' ', start);
    }

    word = str.substr(start, end - start);

    if (word != "") {
        result.push_back(word);
    }

    return result;
}

bool isUnsigned(std::string text) {
    return (text.size() != 0) && (std::find_if(text.begin(), text.end(), [](char a) { return !std::isdigit(a); }) == text.end());
}

/*
* Executes given command if syntax is correct. Stores control flags in bitset flags:
* 1st bit - true if syntax was correct
* 2nd bit - true if program is to be exited
*/
void execute(std::vector<std::string> command, std::bitset<2>& flags) {
    if (command.size() == 0) {
        flags = 0b0;
        return;
    }

    if (command.at(0) == "calcTransfer" || command.at(0) == "t") {
        const std::unordered_set<std::string_view> availableTypes = {
                "square",
                "decoupled",
                "kagome", "kagomePoly", "kagomeAlt"
        };

        if (command.size() < 4 || availableTypes.find(command.at(1)) == availableTypes.end()) {
            flags = 0b0;
            return;
        } //should add a check to make sure the next two arguments are unsigneds

        runCalcTransfer(command.at(1), std::stoul(command.at(2)), std::stoul(command.at(3))); //std::stoui doesnt exist :((
    }
    else if (command.at(0) == "calcEigen" || command.at(0) == "e") {
        const std::unordered_set<std::string_view> availableTypes = {
                "square",
                "decoupled",
                "kagome", "kagomePoly"
        };

        if (command.size() < 3 || availableTypes.find(command.at(1)) == availableTypes.end()) {
            flags = 0b0;
            return;
        } //should add a check to make sure the next two arguments are unsigneds

        if (command.size() > 3) {
            runCalcEigen(command.at(1), std::stoul(command.at(2)), std::stoul(command.at(3))); 
        }
        else {
            runCalcEigen(command.at(1), std::stoul(command.at(2))); //std::stoui doesnt exist :((
        }
    }
    else if (command.at(0) == "polyToFile") {
        if (command.size() < 5) {
            flags = 0b0;
            return;
        }
        for (int i = 1; i < 5; i++) {
            if (!isUnsigned(command.at(i))) {
                flags = 0b0;
                return;
            }
        }

        runPolyToFile(std::stoul(command.at(1)), std::stoul(command.at(2)), std::stoul(command.at(3)), std::stoul(command.at(4)));
    }
    else if (command.at(0) == "exit") {
        flags = 0b11;
        return;
    }
    else {
        flags = 0b0;
        return;
    }

    flags = 0b1;
}

int main() {
    std::bitset<2> executeFlags = 0;

    while (!executeFlags.test(1) && !std::cin.eof()) {
        std::string commandRaw;
        std::vector<std::string> command;

        std::cout << ">";
        std::getline(std::cin, commandRaw);
        command = split(commandRaw);

        execute(command, executeFlags);

        if (!executeFlags.test(0)) {
            std::cout << "Niepoprawna komenda\n";
        }

        std::cout << '\n';
    }
}
