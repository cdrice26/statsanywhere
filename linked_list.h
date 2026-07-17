#ifndef LINKED_LIST_H
#define LINKED_LIST_H
#include <stdlib.h>

typedef struct Node {
    int row;
    int col;
    struct Node *next;
} Node;

Node* create_linked_list();
void push(Node* node, int row, int col);
void free_linked_list(Node* head);

#endif
