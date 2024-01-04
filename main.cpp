#include <iostream>
#include <string.h>

char SkipChars(size_t num) {
  char trash;
  for (size_t i = 0; i < num; ++i) {
    trash = std::getchar();
  }
  return trash;
}

void PrintBack(void** current) {
  if (current == nullptr) {
    std::cout << "error\n";
    return;
  }
  size_t i = 0;
  while (((char*)current[1])[i] != '!') {
    std::cout << ((char*)current[1])[i];
    ++i;
  }
  std::cout << '\n';
}

void** Pop(void** last, size_t& size) {
  if (last == nullptr) {
    return nullptr;
  }
  --size;
  if (last[1] != nullptr) {
    delete[](char*) last[1];
  }
  void** previous = (void**)(last[0]);
  delete[] last;
  return previous;
}

void** Clear(void** current, size_t& size) {
  while (current != nullptr) {
    current = Pop(current, size);
  }
  return current;
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
  string[size] = '!';
  delete[] start;
  return string;
}

void Push(void**& current, size_t& size) {
  SkipChars(3);

  void** new_node = new void* [2];
  new_node[0] = (void*)current;
  new_node[1] = (void*)InputString();
  current = new_node;
  ++size;

  std::cout << "ok\n";
}

void RequestProcessing(void**& current, size_t& size, char symbol) {
  if (symbol == 'b') { //back
    SkipChars(4);
    PrintBack(current);
  }
  if (symbol == 's') { //size
    SkipChars(4);
    std::cout << size << '\n';
  }
  if (symbol == 'c') { //clear
    SkipChars(5);
    current = Clear(current, size);
    std::cout << "ok\n";
  }
  if (symbol == 'p') {
    symbol = getchar();
    if (symbol == 'u') { //push
      Push(current, size);
    } else { //pop
      SkipChars(2);
      PrintBack(current);
      current = Pop(current, size);
    }
  }
}

int main() {
  void** current = nullptr;
  size_t size = 0;

  while (true) {
    char symbol;
    symbol = getchar();
    if (symbol == 'e') { //exit
      SkipChars(4);
      Clear(current, size);
      break;
    }
    RequestProcessing(current, size, symbol);
  }
  std::cout << "bye";
}
