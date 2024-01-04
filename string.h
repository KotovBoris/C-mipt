#include <cstring>
#include <iostream>

char* input_cstring(std::istream& in) {
  size_t capacity = 8;
  size_t size = 0;
  char* start = new char[capacity + 1];

  char symbol = in.get();
  while (std::isspace(symbol) != 0) {
    if (size == capacity) {
      capacity *= 2;
      char* new_start = new char[capacity + 1];
      memmove(new_start, start, size);
      delete[] start;
      start = new_start;
    }

    start[size++] = symbol;

    symbol = in.get();
  }

  start[size] = '\0';
  return start;
}

class String {
public:
  String(const char* cstring) : String(strlen(cstring)) {
    std::copy(cstring, cstring + size_, data_);
  }

  String(size_t length, char symbol) : String(length) {
    std::fill(data_, data_ + size_, symbol);
  }

  String() : String(static_cast<size_t>(0)) {}

  String(const String& other) : String(other.data_) {}

  String& operator=(const String& other) {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    if (capacity_ < size_) {
      delete[] data_;
      data_ = new char[other.capacity_ + 1];
      capacity_ = other.capacity_;
    }
    std::copy(other.data_, other.data_ + size_ + 1, data_); // ÒÐÀÁË Â ÝÒÎÌ, ÕÇ ÏÎ×
    return *this;
  }

  ~String() {
    delete[] data_;
  }

  char& operator[](size_t index) { // out of range => UB
    return data_[index];
  }

  const char& operator[](size_t index) const { // out of range => UB
    return data_[index];
  }

  size_t length() const {
    return size_;
  }

  size_t& size() {
    return size_;
  }

  size_t& capacity() {
    return capacity_;
  }

  void push_back(char symbol) {
    if (size_ == capacity_) {
      capacity_ *= 2;
      this->move_data();
    }
    data_[size_] = symbol;
    ++size_;
    data_[size_] = '\0';
  }

  void pop_back() {
    --size_;
    data_[size_] = '\0';
  }

  char& front() {
    return data_[0];
  }

  char& back() {
    return data_[size_ - 1];
  }

  const char& front() const {
    return data_[0];
  }

  const char& back() const {
    return data_[size_ - 1];
  }

  String& operator+=(char symbol) {
    this->push_back(symbol);
    return *this;
  }

  String& operator+=(const String& second) {
    if (size_ + second.size_ > capacity_) {
      capacity_ = 2 * (size_ + second.size_);
      this->move_data();
    }
    std::copy(second.data_, second.data_ + second.size_ + 1, data_ + size_);
    size_ += second.size_;
    return *this;
  }

  size_t find(const char* substring) const {
    char* csubstring = strstr(data_, substring);
    if (csubstring == nullptr) {
      return size_;
    }
    return csubstring - data_;
  }

  size_t find(const String& substring) const {
    return this->find(substring.data_);
  }

  size_t rfind(const char* substring) const {
    return this->rfind(substring, strlen(substring));
  }

  size_t rfind(const String& substring) const {
    return this->rfind(substring.data_, substring.size_);
  }

  String substr(size_t start, size_t count) const {
    String substr(count);
    std::copy(data_ + start, data_ + start + count, substr.data_);
    return substr;
  }

  bool empty() {
    return (size_ == 0);
  }

  void clear() {
    size_ = 0;
    data_[0] = '\0';
  }

  void shrink_to_fit() {
    capacity_ = size_;
    move_data();
  }

  char*& data() { return data_; }
  char* data() const { return data_; }

private:
  size_t size_;
  size_t capacity_;
  char* data_;

  String(size_t length) : size_(length), capacity_(std::max(size_ * 2, static_cast<size_t>(1))),
    data_(new char[capacity_ + 1]) {
    data_[size_] = '\0';
  }

  void move_data() {
    char* new_data = new char[capacity_ + 1];
    std::copy(data_, data_ + size_ + 1, new_data);
    delete[] data_;
    data_ = new_data;
  }

  size_t rfind(const char* substring, size_t length) const {
    for (long long i = size_ - length; i >= 0; --i) {
      if (memcmp(substring, data_ + i, length) == 0) {
        return i;
      }
    }
    return size_;
  }
};

bool operator==(const String& first, const String& second) {
  return (strcmp(first.data(), second.data()) == 0);
}

bool operator==(const String& first, const char* second) {
  return (strcmp(first.data(), second) == 0);
}

bool operator!=(const String& first, const String& second) {
  return !(first == second);
}

bool operator<(const String& first, const String& second) {
  return (strcmp(first.data(), second.data()) < 0);
}

bool operator>(const String& first, const String& second) {
  return (second < first);
}

bool operator<=(const String& first, const String& second) {
  return !(first > second);
}

bool operator>=(const String& first, const String& second) {
  return !(first < second);
}

String operator+(char symbol, String string) {
  char*& data = string.data();
  size_t& size = string.size();
  size_t& capacity = string.capacity();
  if (size + 1 > capacity) {
    capacity *= 2;
    char* new_data = new char[capacity + 1];
    std::copy(data, data + size + 1, new_data + 1);
    delete[] data;
    data = new_data;
    data[0] = symbol;
    ++size;
    return string;
  }
  memmove(data + 1, data, size + 1);
  data[0] = symbol;
  ++size;
  return string;
}

String operator+(const String& first, const String& second) {
  String copy = first;
  copy += second;
  return copy;
}

String operator+(const String& string, char symbol) {
  String copy = string;
  copy += symbol;
  return copy;
}

std::ostream& operator<<(std::ostream& out, const String& string) {
  //out << string.data();
  for (size_t i = 0; i < string.length(); ++i) {
    out << string.data()[i];
  }
  return out;
}

std::istream& operator>>(std::istream& in, String& string) {
  //char* cstring = input_cstring(in);
  char* cstring = new char[100000];
  in >> cstring;
  string = String(cstring);
  delete[] cstring;
  return in;
}
