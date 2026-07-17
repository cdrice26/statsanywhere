#include "linked_list.h"

Node* create_linked_list() {
    Node *head = (Node*)malloc(sizeof(Node));
    if (head == NULL) return NULL;
    head->next = NULL;
    return head;
}

void push(Node* node, int row, int col) {
    while (node->next != NULL) {
        node = node->next;
    }
    node->row = row;
    node->col = col;
    node->next = (Node*)malloc(sizeof(Node));
    node->next->next = NULL;
}

void free_linked_list(Node* head) {
    Node* temp = head;
    while (temp != NULL) {
        Node *next = temp->next;
        free(temp);
        temp = next;
    }
}
