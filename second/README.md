### Лабораторная работа №2
### «Параллельная реализация решения системы линейных алгебраических уравнений с помощью OpenMP»

#### Замеры времени работы программы при разном количество потоков
_**N = 5 000**_


##### Команда для компиляции. CMakeLists.txt - не рабочий вариант
```bash
gcc -fopenmp main.c sequential_program.c parallel_for_program.c parallel_section_program.c -o main.out -lm`
```
##### Команда для задания количества потоков при исполнении программы
```bash
export OMP_NUM_THREADS=12
./main.out 10000
```