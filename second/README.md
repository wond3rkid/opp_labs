## Лабораторная работа №2
### «Параллельная реализация решения системы линейных алгебраических уравнений с помощью OpenMP»

### some examples
###### Sequential Program
При N = 5000
```c
Time taken: 35.474434 seconds
```
При N = 10000
```c
Time taken for non-parallel: 68.953098 seconds
```


###### Parallel program
При N = 10000
```c
Time taken for parallel: 16.087088 seconds
```

###### Parallel section
При N = 10000
```c
Time taken for 'section' parallel: 52.526029 seconds
```
## 
## 
#### Команда для компиляции. CMakeLists.txt - не рабочий вариант
```bash
 gcc -fopenmp main.c sequential_program.c parallel_for_program.c parallel_section_program.c -o main.out -lm
```