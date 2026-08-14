gcc -g -fsanitize=address -fsanitize=undefined quadrature_weights.c matrices.c vectors.c linked_list.c operations.c special_functions.c main.c -o pfi -lm
