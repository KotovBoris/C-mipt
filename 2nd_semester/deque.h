#include <algorithm>
#include <vector>
#include <iostream>
#include <type_traits>
#include <stdexcept>

template <typename T>
class Deque {
public:
  template <bool IsConst>
  class common_iterator {
  private:
    std::ptrdiff_t ind_chunk = 0;
    std::ptrdiff_t ind_in = 0;
    std::ptrdiff_t chunks_num = 0;
    T** chunks = nullptr;
    T* current_element = nullptr;

    std::ptrdiff_t index() const {
      return ind_chunk * chunk_size + ind_in;
    }

  public:
    using ToReturn = std::conditional_t<IsConst, const T, T>;
    using difference_type = std::ptrdiff_t;
    using value_type = ToReturn;
    using pointer = ToReturn*;
    using reference = ToReturn&;
    using iterator_category = std::random_access_iterator_tag;
    using iterator_type = std::random_access_iterator_tag;

    common_iterator() {};

    common_iterator(Deque& deque) : common_iterator(deque.first) { // return deque.begin()
      if (deque.size() == 0) {
        return;
      }

      chunks_num = deque.chunks_.size();
      chunks = deque.chunks_.data();

      std::ptrdiff_t index = ind_chunk * chunk_size + ind_in;

      if (0 <= index && index < chunks_num * chunk_size) {
        current_element = deque.chunks_[ind_chunk] + ind_in;
      }
      else {
        current_element = nullptr;
      }
    }

    common_iterator& operator+=(std::ptrdiff_t shift) {
      std::ptrdiff_t new_index = index() + shift;

      ind_chunk = new_index / chunk_size;
      ind_in = new_index % chunk_size;

      if (0 <= new_index &&
        new_index < chunk_size * chunks_num) {
        current_element = chunks[ind_chunk] + ind_in;
        return *this;
      }
      current_element = nullptr;
      return *this;
    }

    common_iterator& operator-=(std::ptrdiff_t shift) {
      *this += -shift;
      return *this;
    }

    common_iterator operator+(std::ptrdiff_t shift) const {
      common_iterator copy = *this;
      copy += shift;
      return copy;
    }

    common_iterator operator-(std::ptrdiff_t shift) const {
      common_iterator copy = *this;
      copy -= shift;
      return copy;
    }

    common_iterator& operator++() {
      *this += 1;

      return *this;
    }

    common_iterator operator++(int) {
      common_iterator copy = *this;
      ++(*this);

      return copy;
    }

    common_iterator& operator--() {
      *this -= 1;

      return *this;
    }

    common_iterator operator--(int) {
      common_iterator copy = *this;
      --(*this);

      return copy;
    }

    template<bool Const2>
    std::ptrdiff_t operator-(const common_iterator<Const2>& other) const {
      return index() - other.index();
    }

    bool operator<(const common_iterator& other) const {
      return index() < other.index();
    }

    bool operator<=(const common_iterator& other) const {
      return index() <= other.index();
    }

    bool operator>(const common_iterator& other) const {
      return index() > other.index();
    }

    bool operator>=(const common_iterator& other) const {
      return index() >= other.index();
    }

    bool operator==(const common_iterator& other) const {
      return index() == other.index();
    }

    bool operator!=(const common_iterator& other) const {
      return index() != other.index();
    }

    friend class common_iterator<0>;
    friend class common_iterator<1>;

    operator common_iterator<1>() const {
      common_iterator<1> copy;

      copy.ind_chunk = ind_chunk;
      copy.ind_in = ind_in;
      copy.chunks_num = chunks_num;
      copy.chunks = chunks;
      copy.current_element = current_element;

      return copy;
    }

    reference operator*() const {
      return *current_element;
    }

    pointer operator->() const {
      return current_element;
    }
  };

  using iterator = common_iterator<0>;
  using const_iterator = common_iterator<1>;
  using reverse_iterator = std::reverse_iterator<common_iterator<0>>;
  using const_reverse_iterator = std::reverse_iterator<common_iterator<1>>;

  size_t size() const {
    return size_;
  }

  Deque() : chunks_(0), size_(0) {}

  Deque(std::ptrdiff_t size)
    : chunks_(size / chunk_size + (size % chunk_size != 0)), size_(size) {

    AllocRaw();

    for (auto it = first; it < first + size_; ++it) {
      try {
        new(it.operator->()) T();
      }
      catch (...) {
        Clear(it);
        throw;
      }
    }
  }

  Deque(std::ptrdiff_t size, const T& element)
    : chunks_(size / chunk_size + (size % chunk_size != 0)), size_(size) {

    AllocRaw();

    for (auto it = first; it < first + size_; ++it) {
      try {
        new(it.operator->()) T(element);
      }
      catch (...) {
        Clear(it);
        throw;
      }
    }
  }

  Deque(const Deque& other)
    : chunks_(other.chunks_), first(other.first), size_(other.size_) {

    AllocRaw();

    auto other_it = other.first;
    for (auto it = first; it < first + size_; ++it, ++other_it) {
      try {
        new(&*it) T(*other_it);
      }
      catch (...) {
        Clear(it);
        throw;
      }
    }
  }

  Deque& operator=(const Deque& other) {
    if (this == &other) {
      return *this;
    }

    try {
      Deque copy(other);
      Swap(copy);
    }
    catch (...) {
      throw;
    }

    return *this;
  }

  ~Deque() {
    Clear(end());
  }

  const T& operator[](std::ptrdiff_t index) const {
    return *(first + index);
  }

  T& operator[](std::ptrdiff_t index) {
    return *(first + index);
  }

  const T& at(std::ptrdiff_t index) const {
    if (index < 0 || index >= size_) {
      throw std::out_of_range("out of range in Deque");
    }

    return *(first + index);
  }

  T& at(std::ptrdiff_t index) {
    if (index < 0 || index >= size_) {
      throw std::out_of_range("out of range in Deque");
    }

    return *(first + index);
  }

  void push_back(const T& element) {
    ++size_;

    if ((first + size_ - 1).operator->() == nullptr) {
      std::ptrdiff_t add_chunks = size_ / chunk_size + 1;

      chunks_.reserve(chunks_.size() + add_chunks);
      for (std::ptrdiff_t i = 0; i < add_chunks; ++i) {
        chunks_.push_back(AllocChunk());
      }

      first = iterator(*this);
    }

    try {
      new((first + size_ - 1).operator->()) T(element);
    }
    catch (...) {
      --size_;
      throw;
    }
  }

  void pop_back() {
    --size_;

    (*(first + size_)).~T();
  }

  void push_front(const T& element) {
    ++size_;

    if ((first - 1).operator->() == nullptr) {
      std::ptrdiff_t add_chunks = size_ / chunk_size + 1;

      std::vector<T*> new_chunks(add_chunks);
      for (std::ptrdiff_t i = 0; i < add_chunks; ++i) {
        new_chunks[i] = AllocChunk();
      }

      chunks_.insert(chunks_.begin(), new_chunks.begin(), new_chunks.end());
      first = iterator(*this);
      first += add_chunks * chunk_size;
    }

    --first;

    try {
      new(first.operator->()) T(element);
    }
    catch (...) {
      --size_;
      ++first;
      throw;
    }
  }

  void pop_front() {
    --size_;
    (*first).~T();
    ++first;
  }

  iterator begin() {
    return first;
  }

  const_iterator begin() const {
    return first;
  }

  iterator end() {
    return first + size_;
  }

  const_iterator end() const {
    return first + size_;
  }

  const_iterator cbegin() const {
    return first;
  }

  const_iterator cend() const {
    return first + size_;
  }

  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }

  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }

  reverse_iterator rend() {
    return reverse_iterator(begin());
  }

  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }

  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(end());
  }

  const_reverse_iterator crend() const {
    return const_reverse_iterator(begin());
  }

  void insert(iterator place, const T& element) {
    if (place == end()) {
      push_back(element);
      return;
    }

    T* place_ptr = place.operator->();

    auto last = first + size_ - 1;
    push_back(*last);

    for (auto it = first + size_ - 1; it.operator->() != place_ptr; --it) {
      *it = *(it - 1);
    }

    *place = element;
  }

  void erase(iterator place) {
    auto last = first + size_ - 1;
    for (auto it = place; it != last; ++it) {
      *it = *(it + 1);
    }

    pop_back();
  }

private:
  static const std::ptrdiff_t chunk_size = 2;

  std::vector<T*> chunks_;
  iterator first;
  std::ptrdiff_t size_;

  T* AllocChunk() {
    return reinterpret_cast<T*>(new char[chunk_size * sizeof(T)]);
  }

  void AllocRaw() {
    for (size_t i = 0; i < chunks_.size(); ++i) {
      chunks_[i] = AllocChunk();
    }

    first = iterator(*this);
  }

  void Swap(Deque& other) noexcept {
    chunks_.swap(other.chunks_);
    std::swap(first, other.first);
    std::swap(size_, other.size_);
  }

  void DestructPrefix(iterator end) {
    for (auto it = first; it < end; ++it) {
      (*it).~T();
    }
  }

  void FreeChunks() {
    for (size_t i = 0; i < chunks_.size(); ++i) {
      delete[] reinterpret_cast<char*>(chunks_[i]);
    }
    chunks_.clear();
  }

  void Clear(iterator end) {
    DestructPrefix(end);
    FreeChunks();
  }
};
