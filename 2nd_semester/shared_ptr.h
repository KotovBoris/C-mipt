#include <iostream>
#include <memory>
#include <type_traits>

template<typename T>
class SharedPtr;

template<typename T>
class WeakPtr;

struct BaseControlBlock {
  size_t shared_count = 0;
  size_t weak_count = 0;
  //~BaseControlBlock() = delete;
  virtual void delete_value() = 0;
  virtual void dealloc_cb() = 0;
};

template<typename T, typename Deleter = std::default_delete<T>, typename Alloc = std::allocator<T>>
struct ControlBlockWithDeleter : BaseControlBlock {

  T* ptr;
  Deleter deleter;
  Alloc alloc;

  ControlBlockWithDeleter(T* ptr, Deleter deleter = Deleter(), Alloc alloc = Alloc())
    : ptr(ptr), deleter(std::move(deleter)), alloc(std::move(alloc)) {}

  ~ControlBlockWithDeleter() {
    std::cout << "destructor ControlBlockWithDeleter";
  }

  void delete_value() override {
    if (ptr) {
      deleter(ptr);
    }
    ptr = nullptr;
  }

  using AllocatorTraits = std::allocator_traits<Alloc>;

  using DeleterAllocator = typename AllocatorTraits::template rebind_alloc<Deleter>;
  using DeleterTraits = std::allocator_traits<DeleterAllocator>;

  using AllocAllocator = typename AllocatorTraits::template rebind_alloc<Alloc>;
  using AllocTraits = std::allocator_traits<AllocAllocator>;

  using CBAlloctor = typename AllocatorTraits::template rebind_alloc<ControlBlockWithDeleter>;
  using CBTraits = std::allocator_traits<CBAlloctor>;


  void dealloc_cb() override {
    //delete_value();
    DeleterAllocator deleter_alloc(alloc);
    AllocAllocator alloc_alloc(alloc);
    CBAlloctor cb_alloc(alloc);

    //DeleterTraits::destroy(deleter_alloc, &deleter);
    //AllocTraits::destroy(alloc_alloc, &alloc);

    CBTraits::deallocate(cb_alloc, this, 1);
  }
};

template<typename T, typename Alloc>
struct ControlBlockWithObject : BaseControlBlock {
  T object;
  Alloc alloc;

  template<typename... Args>
  ControlBlockWithObject(Alloc alloc, Args&&... args) : object(std::forward<Args>(args)...), alloc(std::move(alloc)) {}

  ~ControlBlockWithObject() {
    std::cout << "destructor ControlBlockWithObject";
  }

  using AllocatorTraits = std::allocator_traits<Alloc>;

  using AllocAllocator = typename AllocatorTraits::template rebind_alloc<Alloc>;
  using AllocTraits = std::allocator_traits<AllocAllocator>;

  using CBAlloctor = typename AllocatorTraits::template rebind_alloc<ControlBlockWithObject<T, Alloc>>;
  using CBTraits = std::allocator_traits<CBAlloctor>;

  using TAlloctor = typename AllocatorTraits::template rebind_alloc<T>;
  using TTraits = std::allocator_traits<TAlloctor>;

  void delete_value() override {
    TAlloctor t_alloc(alloc);
    TTraits::destroy(t_alloc, &object);
  }

  void dealloc_cb() override {
    //delete_value();
    AllocAllocator alloc_alloc(alloc);
    CBAlloctor cb_alloc(alloc);

    //CBTraits::destroy(cb_alloc, this);     ���������� ��� ����� ����� �� � �� ��� �����

    CBTraits::deallocate(cb_alloc, this, 1);
  }

  //void delete_value() override {
  //  //CBAlloctor cb_alloc(alloc);

  //  //CBTraits::destroy(cb_alloc, this);
  //  //CBTraits::deallocate(cb_alloc, this, 1);
  //}

  //void dealloc_cb() {
  //  CBAlloctor cb_alloc(alloc);

  //  CBTraits::destroy(cb_alloc, this);
  //  CBTraits::deallocate(cb_alloc, this, 1);
  //}
};

template<typename T>
class EnableSharedFromThis {
public:
  SharedPtr<T> shared_from_this();
  WeakPtr<T> weak_from_this();

protected:
  EnableSharedFromThis() = default;
  EnableSharedFromThis(const EnableSharedFromThis&) {}
  EnableSharedFromThis& operator=(const EnableSharedFromThis&) { return *this; }
  ~EnableSharedFromThis() = default;

private:
  WeakPtr<T> weak_ptr_;
};

template<typename T>
class SharedPtr {
public:
  SharedPtr() : control_block(nullptr), ptr(nullptr) {}

  explicit SharedPtr(T* ptr) : control_block(new ControlBlockWithDeleter<T, std::default_delete<T>, std::allocator<T>>(ptr)), ptr(ptr) {
    inc_shared();
    try_enable_shared_from_this(ptr);
  }

  template<typename Y, typename Deleter>
  SharedPtr(Y* ptr, Deleter d) : control_block(create_cb_with_deleter<Y, Deleter, std::allocator<Y>>(ptr, std::move(d))), ptr(ptr) {
    inc_shared();
    try_enable_shared_from_this(ptr);
  }

  template<typename Y, typename Deleter, typename Alloc>
  SharedPtr(Y* ptr, Deleter d, Alloc alloc) : control_block(create_cb_with_deleter<Y, Deleter, Alloc>(ptr, std::move(d), std::move(alloc))), ptr(ptr) {
    inc_shared();
    try_enable_shared_from_this(ptr);
  }

  SharedPtr(const SharedPtr& other) : control_block(other.control_block), ptr(other.ptr) {
    inc_shared();
  }

  SharedPtr(SharedPtr&& other) noexcept : control_block(other.control_block), ptr(other.ptr) {
    other.control_block = nullptr;
    other.ptr = nullptr;
  }

  template<typename Y>
  SharedPtr(const SharedPtr<Y>& other) : control_block(other.control_block), ptr(other.ptr) {
    inc_shared();
  }

  template<typename Y>
  SharedPtr(SharedPtr<Y>&& other) noexcept : control_block(other.control_block), ptr(other.ptr) {
    other.control_block = nullptr;
    other.ptr = nullptr;
  }

  template<typename Y>
  SharedPtr(const SharedPtr<Y>& other, T* ptr) : control_block(other.control_block), ptr(ptr) {
    inc_shared();
  }

  ~SharedPtr() {
    dec_shared();
  }

  SharedPtr& operator=(const SharedPtr& other) {
    if (this != &other) {
      dec_shared();
      control_block = other.control_block;
      ptr = other.ptr;
      inc_shared();
    }
    return *this;
  }

  SharedPtr& operator=(SharedPtr&& other) noexcept {
    if (this != &other) {
      dec_shared();
      control_block = other.control_block;
      ptr = other.ptr;
      other.control_block = nullptr;
      other.ptr = nullptr;
    }
    return *this;
  }

  template<typename Y>
  SharedPtr& operator=(const SharedPtr<Y>& other) {
    dec_shared();
    control_block = other.control_block;
    ptr = other.ptr;
    inc_shared();
    return *this;
  }

  template<typename Y>
  SharedPtr& operator=(SharedPtr<Y>&& other) noexcept {
    dec_shared();
    control_block = other.control_block;
    ptr = other.ptr;
    other.control_block = nullptr;
    other.ptr = nullptr;
    return *this;
  }

  void reset(T* new_ptr = nullptr) {
    //dec_shared();
    //control_block = ptr ? new ControlBlockWithDeleter<T, std::default_delete<T>, std::allocator<T>>(ptr) : nullptr;
    //this->ptr = ptr;
    //if (control_block) {
    //  inc_shared();
    //  try_enable_shared_from_this(ptr);
    //}

    if (new_ptr) {
      *this = SharedPtr(new_ptr);
      return;
    }

    dec_shared();
    control_block = nullptr;
    ptr = nullptr;
  }

  size_t use_count() const {
    return control_block ? control_block->shared_count : 0;
  }

  T* get() const {
    return ptr;
  }

  T& operator*() const {
    return *ptr;
  }

  T* operator->() const {
    return ptr;
  }

  explicit operator bool() const {
    return ptr != nullptr;
  }

  void swap(SharedPtr& other) noexcept {
    std::swap(control_block, other.control_block);
    std::swap(ptr, other.ptr);
  }

private:
  explicit SharedPtr(const WeakPtr<T>& wp) : control_block(wp.control_block), ptr(wp.ptr) {
    inc_shared();
  }

  template<typename Y, typename Deleter, typename Alloc>
  static ControlBlockWithDeleter<Y, Deleter, Alloc>* create_cb_with_deleter(Y* ptr, Deleter d, Alloc alloc = Alloc()) {
    using CB = ControlBlockWithDeleter<Y, Deleter, Alloc>;
    using CBAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<CB>;
    CBAlloc cbAlloc(alloc);
    CB* cb = std::allocator_traits<CBAlloc>::allocate(cbAlloc, 1);

    new (cb) ControlBlockWithDeleter<Y, Deleter, Alloc>(ptr, std::move(d), alloc);

    //cb->deleter = std::move(d);
    //cb->alloc = alloc;
    //cb->shared_count = 0;
    //cb->weak_count = 0;
    //cb->ptr = ptr;

    return cb;
  }

  void inc_shared() {
    if (control_block) {
      ++control_block->shared_count;
    }
  }

  void dec_shared() {
    if (control_block && --control_block->shared_count == 0) {
      
      control_block->delete_value();

      if (control_block->weak_count == 0) {
        control_block->dealloc_cb();
      }

      control_block = nullptr;
      ptr = nullptr;
    }
  }

  template<typename Y>
  void try_enable_shared_from_this(EnableSharedFromThis<Y>* e) {
    if (e) {
      e->weak_ptr_ = *this;
    }
  }

  void try_enable_shared_from_this(...) {}

  template<typename Y>
  friend class SharedPtr;

  template<typename Y>
  friend class WeakPtr;

  BaseControlBlock* control_block;
  T* ptr;

  template<typename U, typename... Args>
  friend SharedPtr<U> makeShared(Args&&... args);

  template<typename U, typename Alloc, typename... Args>
  friend SharedPtr<U> allocateShared(const Alloc& alloc, Args&&... args);
};

template<typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
  using CB = ControlBlockWithObject<T, std::allocator<T>>;
  std::allocator<CB> cbAlloc;
  CB* cb = cbAlloc.allocate(1);
  std::allocator_traits<std::allocator<CB>>::construct(cbAlloc, cb, std::allocator<T>(), std::forward<Args>(args)...);
  SharedPtr<T> sp;
  sp.control_block = cb;
  sp.ptr = &cb->object;
  sp.inc_shared();
  sp.try_enable_shared_from_this(sp.ptr);
  return sp;
}

template<typename T, typename Alloc, typename... Args>
SharedPtr<T> allocateShared(const Alloc& alloc, Args&&... args) {
  using CB = ControlBlockWithObject<T, Alloc>;
  using CBAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<CB>;
  CBAlloc cbAlloc(alloc);
  CB* cb = std::allocator_traits<CBAlloc>::allocate(cbAlloc, 1);
  std::allocator_traits<CBAlloc>::construct(cbAlloc, cb, alloc, std::forward<Args>(args)...);
  SharedPtr<T> sp;
  sp.control_block = cb;
  sp.ptr = &cb->object;
  sp.inc_shared();
  sp.try_enable_shared_from_this(sp.ptr);
  return sp;
}

template<typename T>
class WeakPtr {
public:
  WeakPtr() : control_block(nullptr), ptr(nullptr) {}

  WeakPtr(const SharedPtr<T>& sp) : control_block(sp.control_block), ptr(sp.ptr) {
    inc_weak();
  }

  template<typename Y>
  WeakPtr(const SharedPtr<Y>& sp) : control_block(sp.control_block), ptr(sp.ptr) {
    inc_weak();
  }

  WeakPtr(const WeakPtr& other) : control_block(other.control_block), ptr(other.ptr) {
    inc_weak();
  }

  template<typename Y>
  WeakPtr(const WeakPtr<Y>& other) : control_block(other.control_block), ptr(other.ptr) {
    inc_weak();
  }

  template<typename Y>
  WeakPtr(WeakPtr<Y>&& other) noexcept : control_block(other.control_block), ptr(other.ptr) {
    other.control_block = nullptr;
    other.ptr = nullptr;
  }

  ~WeakPtr() {
    dec_weak();
  }

  WeakPtr& operator=(const WeakPtr& other) {
    if (this != &other) {
      dec_weak();
      control_block = other.control_block;
      ptr = other.ptr;
      inc_weak();
    }
    return *this;
  }

  template<typename Y>
  WeakPtr& operator=(const WeakPtr<Y>& other) {
    dec_weak();
    control_block = other.control_block;
    ptr = other.ptr;
    inc_weak();
    return *this;
  }

  WeakPtr& operator=(WeakPtr&& other) noexcept {
    if (this != &other) {
      dec_weak();
      control_block = other.control_block;
      ptr = other.ptr;
      other.control_block = nullptr;
      other.ptr = nullptr;
    }
    return *this;
  }

  template<typename Y>
  WeakPtr& operator=(WeakPtr<Y>&& other) noexcept {
    dec_weak();
    control_block = other.control_block;
    ptr = other.ptr;
    other.control_block = nullptr;
    other.ptr = nullptr;
    return *this;
  }

  bool expired() const {
    return !control_block || control_block->shared_count == 0;
  }

  SharedPtr<T> lock() const {
    return expired() ? SharedPtr<T>() : SharedPtr<T>(*this);
  }

  size_t use_count() const {
    return control_block ? control_block->shared_count : 0;
  }

private:
  void inc_weak() {
    if (control_block) {
      ++control_block->weak_count;
    }
  }

  void dec_weak() {
    if (control_block && --control_block->weak_count == 0 && control_block->shared_count == 0) {
      
      control_block->dealloc_cb();

      control_block = nullptr;
      ptr = nullptr;
    }
  }

  template<typename>
  friend class SharedPtr;

  template<typename>
  friend class WeakPtr;

  BaseControlBlock* control_block;
  T* ptr;
};

template<typename T>
SharedPtr<T> EnableSharedFromThis<T>::shared_from_this() {
  return SharedPtr<T>(weak_from_this());
}

template<typename T>
WeakPtr<T> EnableSharedFromThis<T>::weak_from_this() {
  return weak_ptr_;
}
