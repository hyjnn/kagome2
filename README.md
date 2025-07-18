Na linuxie używałem polecenia
```
g++ -std=c++20 source/*.cpp -O3 -I"/path/to/eigen" -lquadmath -o kagome
```
żeby skompilować kod (należy je wykonać będąc w głównym katalogu repozytorium).
Ścieżkę "/path/to/eigen" należy zastąpić ścieżką do folderu zawierającego bibiliotekę Eigen (powinien w nim znajdować się folder Eigen, w którym zawarte są nagłówki Eigena)
Program domyślnie wykorzystuje zmienne o podwójnej precyzji przy niektórych obliczeniach, aby pozwolić na wykorzystanie biblioteki Intel MKL do przyspieszenia obliczeń. Kompilując program z dodatkową flagą USE\_FLOAT128\_ALL (tj. dodając opcję -DUSE\_FLOAT128\_ALL) możemy przełączyć program do działania z poczwórnej precyzją we wszystkich obliczeniach.
