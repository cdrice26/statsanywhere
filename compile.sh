gcc -g -fsanitize=address -fsanitize=undefined quadrature_weights.c normal.c matrices.c vectors.c linked_list.c -o pfi -lm
