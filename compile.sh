gcc -g -fsanitize=address -fsanitize=undefined quadrature_weights.c main.c matrices.c vectors.c linked_list.c operations.c -o pfi -lm
