### Лабораторная работа №1
### «Параллельная реализация решения системы линейных алгебраических уравнений с помощью MPI»

Исполнение файла с size = 4
```bash
mpiexec -np 4 ./first
```
Исполнение файла под valgrind для определения ошибок
```bash
mpiexec -np 4 valgrind ./first
```