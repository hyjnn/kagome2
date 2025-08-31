## Użycie

Program jest aplikacją działającą z poziomu konsoli. Po wywołaniu, program oczekuje na podanie przez użytkownika jednej z następujących komend:
- `calcTransfer type n m`, `t type n m` - Liczy entropię resztkową sieci type o wymiarze n x m (obsługiwane wartości type to: square, decoupled, kagome, kagomePoly) licząc ślad macierzy transferu. Dla mniejszych sieci wyświetla macierz transferu.
- `calcEigen type n iter`, `e type n iter` - Liczy entropię resztkową sieci type o wymiarze n x inf (obsługowiane wartości type jak wyżej). Argument iter określa liczbę iteracji metody potęgowej stosowanej do obliczenia maksymalnej wartości własnej.
- `exit` - Zamyka program.

## Kompilacja

Na linuxie (kompilator gcc) używałem polecenia
```
g++ -std=c++20 source/*.cpp -O3 -I"/path/to/eigen" -lquadmath -o kagome
```
żeby skompilować kod (należy je wykonać będąc w głównym katalogu repozytorium). Program można również skompilować kompilatorem MSVC analogiczną komendą (bez -lquadmath). Program nie wspiera wtedy obliczeń w poczwórnej precyzji.
Ścieżkę "/path/to/eigen" należy zastąpić ścieżką do folderu zawierającego bibiliotekę Eigen (powinien w nim znajdować się folder Eigen, w którym zawarte są nagłówki Eigena)
Program domyślnie wykorzystuje zmienne o podwójnej precyzji przy niektórych obliczeniach. Aby program wykonywał wszystkie obliczenia w poczwórnej precyzji, należy skompilować z flagą USE\_FLOAT128\_ALL:
```
g++ -std=c++20 source/*.cpp -O3 -I"/path/to/eigen" -lquadmath -o kagome -DUSE_FLOAT128_ALL
```
Bez użycia flagi program można dodatkowo skompilować z użyciem biblioteki Intel MKL, co skraca czas obliczeń w podwójnej precyzji.
