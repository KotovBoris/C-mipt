#include <memory>
#include <type_traits>

template<size_t N>
class StackStorage {
public:
  StackStorage() = default;
  StackStorage(const StackStorage&) = delete;
  StackStorage& operator=(const StackStorage&) = delete;

  char& operator[](size_t index) {
    return data_[index];
  }

  template<typename T>
  T* reserve(size_t count) {
    std::align(alignof(T), sizeof(T) * count, top_, free_space_);
    T* reserved = reinterpret_cast<T*>(top_);

    top_ = reinterpret_cast<void*>(reserved + count);
    free_space_ -= sizeof(T) * count;

    return reserved;
  }

private:
  char data_[N];
  void* top_ = reinterpret_cast<void*>(data_);
  size_t free_space_ = N;
};

template<typename T, size_t N>
class StackAllocator {
public:
  using value_type = T;

  StackAllocator() = default;
  explicit StackAllocator(StackStorage<N>& storage) : storage_(&storage) {}

  template<typename U, size_t M>
  friend class StackAllocator;

  template<typename U>
  StackAllocator(const StackAllocator<U, N>& other) : storage_(other.storage_) {}

  template<typename U>
  StackAllocator& operator=(const StackAllocator<U, N>& other) {
    storage_ = other.storage_;
    return *this;
  }

  T* allocate(size_t count) {
    return storage_->template reserve<T>(count);
  }

  void deallocate(T*, size_t) {}

  template<typename U, size_t M>
  bool operator==(const StackAllocator<U, M>& other) const {
    return storage_ == other.storage_;
  }

  template<typename U, size_t M>
  bool operator!=(const StackAllocator<U, M>& other) const {
    return storage_ != other.storage_;
  }

  template <typename U>
  struct rebind {
    using other = StackAllocator<U, N>;
  };

private:
  StackStorage<N>* storage_ = nullptr;
};

template<typename T, typename Allocator = std::allocator<T>>
class List {
private:
  struct BaseNode {
    BaseNode* previous_;
    BaseNode* next_;
  };

  struct Node : BaseNode {
    Node(BaseNode* previous, BaseNode* next) : BaseNode{ previous, next } {}
    Node(BaseNode* previous, BaseNode* next, const T& value)
      : BaseNode{ previous, next }, value_(value) {}

    T value_;
  };

  using AllocTraits = std::allocator_traits<Allocator>;
  using NodeAllocator = typename AllocTraits::template rebind_alloc<Node>;
  using NodeTraits = std::allocator_traits<NodeAllocator>;

public:
  template <bool IsConst>
  class common_iterator {
  public:
    friend class List;

    using ToReturn = std::conditional_t<IsConst, const T, T>;
    using difference_type = std::ptrdiff_t;
    using value_type = ToReturn;
    using pointer = ToReturn*;
    using reference = ToReturn&;
    using iterator_category = std::bidirectional_iterator_tag;

    common_iterator() = default;
    explicit common_iterator(const List& list) : ptr_(list.fake_.next_) {}

    common_iterator& operator++() {
      ptr_ = ptr_->next_;
      return *this;
    }

    common_iterator operator++(int) {
      common_iterator copy = *this;
      ++(*this);
      return copy;
    }

    common_iterator& operator--() {
      ptr_ = ptr_->previous_;
      return *this;
    }

    common_iterator operator--(int) {
      common_iterator copy = *this;
      --(*this);
      return copy;
    }

    bool operator==(const common_iterator& other) const {
      return ptr_ == other.ptr_;
    }

    bool operator!=(const common_iterator& other) const {
      return ptr_ != other.ptr_;
    }

    friend class common_iterator<false>;
    friend class common_iterator<true>;

    operator common_iterator<true>() const {
      common_iterator<true> copy;
      copy.ptr_ = ptr_;
      return copy;
    }

    reference operator*() const {
      Node* real_ptr = reinterpret_cast<Node*>(ptr_);
      return real_ptr->value_;
    }

    pointer operator->() const {
      return &*(*this);
    }

  private:
    BaseNode* ptr_;
  };

  using iterator = common_iterator<false>;
  using const_iterator = common_iterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  List() = default;
  explicit List(size_t size) : allocator_() { CreateDefault(size); }
  List(size_t size, const T& element) : allocator_() { CreateSpecific(size, element); }
  List(const Allocator& allocator) : allocator_(allocator) {}
  List(size_t size, const Allocator& allocator) : allocator_(allocator) { CreateDefault(size); }
  List(size_t size, const T& element, const Allocator& allocator) : allocator_(allocator) { CreateSpecific(size, element); }
  List(const List& other) : allocator_(AllocTraits::select_on_container_copy_construction(other.get_allocator())) { Fill(other); }

  ~List() { clear(); }

  List& operator=(const List& other) {
    bool copy_alloc = AllocTraits::propagate_on_container_copy_assignment::value;
    Allocator alloc = copy_alloc ? other.get_allocator() : get_allocator();
    List copy(alloc);
    copy.Fill(other);
    Swap(copy);
    return *this;
  }

  Allocator get_allocator() const {
    return allocator_;
  }

  size_t size() const {
    return size_;
  }

  void clear() {
    while (size() > 0) {
      pop_back();
    }
  }

  iterator begin() { return iterator(*this); }
  const_iterator begin() const { return iterator(*this); }
  iterator end() { iterator fake_node = begin(); --fake_node; return fake_node; }
  const_iterator end() const { const_iterator fake_node = cbegin(); --fake_node; return fake_node; }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }
  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator crend() const { return const_reverse_iterator(begin()); }

  void push_front(const T& element) { insert(begin(), element); }
  void push_back(const T& element) {
    insert(end(), element);
  }

  void pop_front() {
    erase(begin());
  }

  void pop_back() {
    erase(last());
  }

  void insert(const_iterator it, const T& element) {
    Node* new_node = AllocNode();
    NodeAllocator alloc(get_allocator());
    try {
      NodeTraits::construct(alloc, new_node, it.ptr_->previous_, it.ptr_, element);
    }
    catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }

    it.ptr_->previous_->next_ = new_node;
    it.ptr_->previous_ = new_node;

    ++size_;
  }

  void erase(const_iterator it) {
    it.ptr_->previous_->next_ = it.ptr_->next_;
    it.ptr_->next_->previous_ = it.ptr_->previous_;
    NodeAllocator alloc(get_allocator());
    NodeTraits::destroy(alloc, reinterpret_cast<Node*>(it.ptr_));
    NodeTraits::deallocate(alloc, reinterpret_cast<Node*>(it.ptr_), 1);

    --size_;
  }

private:
  Node* AllocNode() {
    NodeAllocator alloc(get_allocator());
    Node* ptr = NodeTraits::allocate(alloc, 1);
    return ptr;
  }

  void PushDefault() {
    iterator it = end();
    Node* new_node = AllocNode();
    NodeAllocator alloc(get_allocator());
    try {
      NodeTraits::construct(alloc, new_node, it.ptr_->previous_, it.ptr_);
    }
    catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }

    it.ptr_->previous_->next_ = new_node;
    it.ptr_->previous_ = new_node;

    ++size_;
  }

  void CreateDefault(size_t number) {
    try {
      for (size_t i = 0; i < number; ++i) {
        PushDefault();
      }
    }
    catch (...) {
      clear();
      throw;
    }
  }

  void CreateSpecific(size_t number, const T& element) {
    try {
      for (size_t i = 0; i < number; ++i) {
        push_back(element);
      }
    }
    catch (...) {
      clear();
      throw;
    }
  }

  void Swap(List& other) {
    std::swap(size_, other.size_);
    std::swap(fake_.previous_, other.fake_.previous_);
    std::swap(fake_.next_->previous_, other.fake_.next_->previous_);
    std::swap(fake_.next_, other.fake_.next_);
    std::swap(fake_.previous_->next_, other.fake_.previous_->next_);

    Allocator tmp = get_allocator();
    allocator_ = other.get_allocator();
    other.allocator_ = tmp;
  }

  iterator last() {
    iterator before_end = end();
    --before_end;
    return before_end;
  }

  void Fill(const List& other) {
    try {
      for (const auto& element : other) {
        push_back(element);
      }
    }
    catch (...) {
      clear();
      throw;
    }
  }

  List(Allocator alloc, const List& other) : Allocator(alloc) {
    Fill(other);
  }

  BaseNode fake_ = { &fake_, &fake_ };
  size_t size_ = 0;
  [[no_unique_address]] Allocator allocator_;
};
