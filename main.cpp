#include <iostream>
#include <string.h>

struct Node {
  char* string;
  Node* previous = nullptr;
};

void SkipChars(size_t num) {
  for (size_t i = 0; i < num; ++i) {
    getchar();
  }
}

void PrintBack(Node* current) {
  if (current == nullptr) {
    std::cout << "error\n";
    return;
  }

  std::cout << current->string << '\n';
}

Node* Pop(Node* last, size_t& size) {
  if (last == nullptr) {
    return nullptr;
  }

  --size;

  Node* previous = last->previous;
  delete[] last->string;
  delete[] last;

  return previous;
}

void Clear(Node* current, size_t& size) {
  while (current != nullptr) {
    current = Pop(current, size);
  }
}

char* InputString() {
  size_t capacity = 8;
  size_t size = 0;
  char* start = new char[capacity];

  char symbol = getchar();
  while (symbol != '\n') {
    if (size == capacity) {
      capacity *= 2;
      char* new_start = new char[capacity];
      memmove(new_start, start, size);
      delete[] start;
      start = new_start;
    }

    start[size++] = symbol;

    symbol = getchar();
  }

  char* string = new char[size + 1];
  memmove(string, start, size);
  string[size] = '\0';
  delete[] start;
  return string;
}

void Push(Node*& current, size_t& size) {
  Node* new_node = new Node{InputString(), current};
  current = new_node;
  ++size;
}

void RequestProcessing(Node*& current, size_t& size, char symbol) {
  if (symbol == 'b') {  // back
    SkipChars(4);
    PrintBack(current);
  }
  if (symbol == 's') {  // size
    SkipChars(4);
    std::cout << size << '\n';
  }
  if (symbol == 'c') {  // clear
    SkipChars(5);
    Clear(current, size);
    current = nullptr;
    std::cout << "ok\n";
  }
  if (symbol == 'p') {
    symbol = getchar();
    if (symbol == 'u') {  // push
      SkipChars(3);
      Push(current, size);
      std::cout << "ok\n";
    } else {  // pop
      SkipChars(2);
      PrintBack(current);
      current = Pop(current, size);
    }
  }
}

int main() {
  Node* current = nullptr;
  size_t size = 0;

  while (true) {
    char symbol;
    symbol = getchar();
    if (symbol == 'e') {  // exit
      SkipChars(4);
      Clear(current, size);
      break;
    }
    RequestProcessing(current, size, symbol);
  }
  std::cout << "bye";
}
