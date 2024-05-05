#include <cstring>
#include <iostream>

class String {
 public:
  String(const char* cstring) : String(strlen(cstring)) {
    std::copy(cstring, cstring + size_, data_);
  }

  String(char symbol) : String(static_cast<size_t>(1)) { data_[0] = symbol; }

  String(size_t length, char symbol) : String(length) {
    std::fill(data_, data_ + size_, symbol);
  }

  String() : String(static_cast<size_t>(0)) {}

  String(const String& other) : String(other.data_) {}

  String& operator=(const String& other);

  ~String() { delete[] data_; }

  char& operator[](size_t index) { return data_[index]; }

  const char& operator[](size_t index) const { return data_[index]; }

  size_t length() const { return size_; }

  size_t size() const { return size_; }

  size_t capacity() const { return capacity_; }

  void push_back(char symbol);

  void pop_back();

  char& front() { return data_[0]; }

  char& back() { return data_[size_ - 1]; }

  const char& front() const { return data_[0]; }

  const char& back() const { return data_[size_ - 1]; }

  String& operator+=(const String& second);

  size_t find(const String& substring) const;

  size_t rfind(const String& substring) const;

  String substr(size_t start, size_t count) const;

  bool empty() const { return (size_ == 0); }

  void clear();

  void shrink_to_fit() { reserve(size_); }

  char*& data() { return data_; }
  char* data() const { return data_; }

 private:
  size_t size_;
  size_t capacity_;
  char* data_;

  String(size_t length);

  void reserve(size_t capacity);
};

String::String(size_t length)
    : size_(length),
      capacity_(std::max(size_ * 2, static_cast<size_t>(1))),
      data_(new char[capacity_ + 1]) {
  data_[size_] = '\0';
}

void String::reserve(size_t capacity) {
  capacity_ = capacity;
  char* new_data = new char[capacity_ + 1];
  std::copy(data_, data_ + size_ + 1, new_data);
  delete[] data_;
  data_ = new_data;
}

String& String::operator=(const String& other) {
  if (this == &other) {
    return *this;
  }
  size_ = other.size_;
  if (capacity_ < size_) {
    reserve(other.capacity_);
  }
  std::copy(other.data_, other.data_ + size_ + 1, data_);
  return *this;
}

void String::push_back(char symbol) {
  if (size_ == capacity_) {
    reserve(capacity_ * 2);
  }
  data_[size_] = symbol;
  ++size_;
  data_[size_] = '\0';
}

void String::pop_back() {
  --size_;
  data_[size_] = '\0';
}

String& String::operator+=(const String& second) {
  size_t new_size = size_ + second.size_;
  if (new_size > capacity_) {
    reserve(2 * new_size);
  }

  std::copy(second.data_, second.data_ + second.size_ + 1, data_ + size_);
  std::copy(second.data_, second.data_ + second.size_ + 1, data_ + size_);
  size_ = new_size;
  return *this;
}

size_t String::find(const String& substring) const {
  char* csubstring = strstr(data_, substring.data_);
  if (csubstring == nullptr) {
    return size_;
  }
  return csubstring - data_;
}

size_t String::rfind(const String& substring) const {
  for (long long i = size_ - substring.size_; i >= 0; --i) {
    if (memcmp(substring.data_, data_ + i, substring.size_) == 0) {
      return i;
    }
  }
  return size_;
}

String String::substr(size_t start, size_t count) const {
  String substr(count);
  std::copy(data_ + start, data_ + start + count, substr.data_);
  return substr;
}

void String::clear() {
  size_ = 0;
  data_[0] = '\0';
}

bool operator==(const String& first, const String& second) {
  return (strcmp(first.data(), second.data()) == 0);
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

String operator+(const String& first, const String& second) {
  String copy = first;
  copy += second;
  return copy;
}

std::ostream& operator<<(std::ostream& out, const String& string) {
  string.data()[string.size()] = '\0';
  out << string.data();
  return out;
}

std::istream& operator>>(std::istream& in, String& string) {
  char symbol;
  while (!in.eof() && in.get(symbol) && !std::isspace(symbol)) {
    string.push_back(symbol);
  }
  return in;
}
