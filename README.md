Na linuxie używałem polecenia
```
g++ -std=c++20 source/*.cpp -O3 -I"/path/to/eigen" -lquadmath -o kagome
```
żeby skompilować kod (należy je wykonać będąc w głównym katalogu repozytorium).
Ścieżkę "/path/to/eigen" należy zastąpić ścieżką do folderu zawierającego bibiliotekę Eigen (powinien w nim znajdować się folder Eigen, w którym zawarte są nagłówki Eigena)
