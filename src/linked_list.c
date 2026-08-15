// linked_list.c
#include "linked_list.h"

/**
 * @brief Creates an empty linked list head node.
 *
 * Allocates and returns a head Node whose data pointer is NULL and next
 * pointer is NULL. Caller is responsible for freeing the list with
 * free_linked_list when no longer needed.
 *
 * @return Pointer to newly-allocated head Node, or NULL on allocation failure
 */
Node* create_linked_list(void) {
    Node *head = malloc(sizeof(Node));
    if (head == NULL) return NULL;
    head->data = NULL;
    head->next = NULL;
    return head;
}

/**
 * @brief Appends data to the end of the linked list starting at node.
 *
 * Traverses from the provided node to the tail, stores data into the
 * current tail node and allocates a new empty tail node. If node is NULL
 * or allocation fails, the function returns 0 to indicate failure.
 *
 * @param node Starting node (typically the head)
 * @param data Pointer to data to store in the new element
 * @return 1 on success, 0 on failure
 */
int push(Node *node, void *data) {
    if (node == NULL) return 0;
    while (node->next != NULL) {
        node = node->next;
    }
    node->data = data;
    node->next = malloc(sizeof(Node));
    if (node->next == NULL) return 0;
    node->next->data = NULL;
    node->next->next = NULL;
    return 1;
}

/**
 * @brief Frees a linked list and optionally its element data.
 *
 * Walks the list starting at head, calling free_data on each non-NULL
 * data pointer if free_data is provided, then frees each Node. Passing
 * NULL for free_data skips freeing element data (useful for stack-allocated
 * or externally-managed data).
 *
 * @param head Pointer to the head of the list (may be NULL)
 * @param free_data Optional callback to free each element's data; pass
 *                  NULL if element data does not require freeing
 */
void free_linked_list(Node *head, void (*free_data)(void*)) {
    Node *temp = head;
    while (temp != NULL) {
        Node *next = temp->next;
        if (free_data && temp->data) free_data(temp->data);
        free(temp);
        temp = next;
    }
}

/**
 * @brief Returns the number of stored elements in the list.
 *
 * Counts nodes with non-empty data fields, treating the final empty tail
 * node (with next == NULL and data == NULL) as not containing an element.
 *
 * @param head Pointer to the head of the list
 * @return Number of elements stored in the list
 */
int length(Node *head) {
    int len = 0;
    Node *temp = head;
    while (temp->next != NULL) {
        len++;
        temp = temp->next;
    }
    return len;
}
