#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>

template<typename T, typename Allocator = std::allocator<T>>
class List : Allocator {
private:
  using AllocTraits = std::allocator_traits<Allocator>;
  using ValueAllocator = typename AllocTraits::template rebind_alloc<T>;
  using ValueTraits = std::allocator_traits<ValueAllocator>;

  struct BaseNode {
    BaseNode* previous;
    BaseNode* next;
  };

  struct Node : BaseNode {
    List* dad;
    T* value_ptr;

    template<typename... Args>
    Node(BaseNode* previous, BaseNode* next, List* dad, Args&&... args)
      : BaseNode{ previous, next }, dad(dad) {
      ValueAllocator alloc(dad->get_allocator());

      value_ptr = ValueTraits::allocate(alloc, 1);
      ValueTraits::construct(alloc, value_ptr, std::forward<Args>(args)...);
    }

    ~Node() {
      ValueAllocator alloc(dad->get_allocator());

      ValueTraits::destroy(alloc, value_ptr);
      ValueTraits::deallocate(alloc, value_ptr, 1);
    }
  };

  BaseNode fake = { &fake, &fake };
  size_t size_ = 0;
  using NodeAllocator = typename AllocTraits::template rebind_alloc<Node>;
  using NodeTraits = std::allocator_traits<NodeAllocator>;

public:
  template <bool IsConst>
  class common_iterator {
  private:
    BaseNode* ptr;

  public:
    friend class List;

    using ToReturn = std::conditional_t<IsConst, const T, T>;
    using difference_type = std::ptrdiff_t;
    using value_type = ToReturn;
    using pointer = ToReturn*;
    using reference = ToReturn&;
    using iterator_category = std::bidirectional_iterator_tag;

    common_iterator() {};

    common_iterator(const List& list) : ptr(list.fake.next) {} // return list.begin()

    common_iterator& operator++() {
      ptr = ptr->next;

      return *this;
    }

    common_iterator operator++(int) {
      common_iterator copy = *this;
      ++(*this);

      return copy;
    }

    common_iterator& operator--() {
      ptr = ptr->previous;

      return *this;
    }

    common_iterator operator--(int) {
      common_iterator copy = *this;
      --(*this);

      return copy;
    }

    bool operator==(const common_iterator& other) const {
      return ptr == other.ptr;
    }

    friend class common_iterator<0>;
    friend class common_iterator<1>;

    operator common_iterator<1>() const {
      common_iterator<1> copy;
      copy.ptr = ptr;

      return copy;
    }

    explicit operator common_iterator<0>() const {
      common_iterator<0> copy;
      copy.ptr = ptr;

      return copy;
    }

    reference operator*() const {
      Node* real_ptr = reinterpret_cast<Node*>(ptr);
      return *(real_ptr->value_ptr);
    }

    pointer operator->() const {
      Node* real_ptr = reinterpret_cast<Node*>(ptr);
      return real_ptr->value_ptr;
    }
  };

  using iterator = common_iterator<0>;
  using const_iterator = common_iterator<1>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() {
    return iterator(*this);
  }

  const_iterator begin() const {
    return iterator(*this);
  }

  iterator end() {
    iterator fake_node = begin();
    --fake_node;
    return fake_node;
  }

  const_iterator end() const {
    const_iterator fake_node = cbegin();
    --fake_node;
    return fake_node;
  }

  const_iterator cbegin() const {
    return begin();
  }

  const_iterator cend() const {
    return end();
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

  List() {}

  Allocator get_allocator() const {
    return static_cast<Allocator>(*this);
  }

  List(size_t size) {
    CreateDefault(size);
  }

  List(size_t size, const T& element) {
    CreateSpecific(size, element);
  }

  List(const Allocator& allocator) : Allocator(allocator) {}

  List(size_t size, const Allocator& allocator) : Allocator(allocator) {
    CreateDefault(size);
  }

  List(size_t size, const T& element, const Allocator& allocator)
    : Allocator(allocator) {
    CreateSpecific(size, element);
  }

  List(const List& other)
    : Allocator(AllocTraits::select_on_container_copy_construction(other.get_allocator())) {
    Fill(other);
  }

  List(List&& other)
    : Allocator(other.get_allocator()) {
    MoveNodes(std::move(other));
  }

  void clear() {
    while (size() > 0) {
      pop_back();
    }
  }

  ~List() {
    clear();
  }

  List& operator=(const List& other) {
    bool copy_alloc = AllocTraits::propagate_on_container_copy_assignment::value;
    Allocator alloc = copy_alloc ? other.get_allocator() : get_allocator();
    List copy(alloc, other);
    swap(copy);
    return *this;
  }

  List& operator=(List&& other) {
    if (this == &other) {
      return *this;
    }

    clear();

    bool move_alloc = AllocTraits::propagate_on_container_move_assignment::value;
    if (move_alloc) {
      Allocator& this_alloc = static_cast<Allocator&>(*this);
      this_alloc = std::move(other.get_allocator());
    }

    MoveNodes(std::move(other));

    return *this;
  }

  size_t size() const {
    return size_;
  }

  void push_front(const T& element) {
    insert(begin(), element);
  }

  void push_back(const T& element) {
    insert(end(), element);
  }

  void push_front(T&& element) {
    insert(begin(), std::move(element));
  }

  void push_back(T&& element) {
    insert(end(), std::move(element));
  }

  void pop_front() {
    erase(begin());
  }

  void pop_back() {
    erase(last());
  }

  iterator insert(const_iterator it, const T& element) {
    Node* new_node = AllocNode();
    NodeAllocator alloc(get_allocator());
    NodeTraits::construct(alloc, new_node, it.ptr->previous, it.ptr,
      this, element);

    it.ptr->previous->next = new_node;
    it.ptr->previous = new_node;

    ++size_;

    return iterator(--it);
  }

  iterator insert(const_iterator it, T&& element) {
    Node* new_node = AllocNode();
    NodeAllocator alloc(get_allocator());
    NodeTraits::construct(alloc, new_node, it.ptr->previous, it.ptr,
      this, std::move(element));

    it.ptr->previous->next = new_node;
    it.ptr->previous = new_node;

    ++size_;

    return iterator(--it);
  }

  template<typename... Args>
  iterator emplace(const_iterator pos, Args&&... args) {
    Node* new_node = AllocNode();
    try {
      NodeAllocator alloc(get_allocator());
      NodeTraits::construct(alloc, new_node, pos.ptr->previous, pos.ptr,
        this, std::forward<Args>(args)...);
    }
    catch (...) {
      NodeAllocator alloc(get_allocator());
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }

    new_node->previous->next = new_node;
    new_node->next->previous = new_node;

    ++size_;

    return iterator(--pos);
  }

  iterator erase(const_iterator it) {
    iterator next = iterator(it);
    ++next;

    it.ptr->previous->next = it.ptr->next;
    it.ptr->next->previous = it.ptr->previous;
    NodeAllocator alloc(get_allocator());
    NodeTraits::destroy(alloc, reinterpret_cast<Node*>(it.ptr));
    NodeTraits::deallocate(alloc, reinterpret_cast<Node*>(it.ptr), 1);

    --size_;

    return next;
  }

  void swap(List& other) {
    std::swap(size_, other.size_);
    std::swap(fake.previous, other.fake.previous);
    std::swap(fake.next->previous, other.fake.next->previous);
    std::swap(fake.next, other.fake.next);
    std::swap(fake.previous->next, other.fake.previous->next);

    Allocator tmp = get_allocator();
    Allocator& this_casted = static_cast<Allocator&>(*this);
    this_casted = other.get_allocator();
    Allocator& other_casted = static_cast<Allocator&>(other);
    other_casted = tmp;
  }

  void MoveNode(iterator node, iterator place) {
    if (node == place || node == end()) {
      return;
    }

    node.ptr->previous->next = node.ptr->next;
    node.ptr->next->previous = node.ptr->previous;

    node.ptr->previous = place.ptr->previous;
    node.ptr->next = place.ptr;
    place.ptr->previous->next = node.ptr;
    place.ptr->previous = node.ptr;
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
    NodeTraits::construct(alloc, new_node, it.ptr->previous, it.ptr, this);

    it.ptr->previous->next = new_node;
    it.ptr->previous = new_node;

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

  void MoveNodes(List&& other) {
    fake.next = other.fake.next;
    fake.previous = other.fake.previous;
    other.fake.next->previous = &fake;
    other.fake.previous->next = &fake;
    size_ = other.size_;

    other.fake.next = &other.fake;
    other.fake.previous = &other.fake;
    other.size_ = 0;
  }
};

template<typename Key, typename Value, typename Hash = std::hash<Key>,
  typename Equal = std::equal_to<Key>, typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap : private Hash, private Equal, private Alloc {
public:
  using NodeType = std::pair<Key, Value>;

private:
  using AllocTraits = std::allocator_traits<Alloc>;
  using ListType = List<NodeType, Alloc>;

  size_t size_;
  ListType list_;
  double load_factor_;
  double max_load_factor_;
  std::vector<typename ListType::iterator> buckets_;

public:
  using iterator = typename ListType::iterator;
  using const_iterator = typename ListType::const_iterator;

  UnorderedMap() : size_(0), load_factor_(0.0), max_load_factor_(1.0) {
    buckets_.resize(16, list_.end());
  }

  UnorderedMap(const UnorderedMap& other)
    : Hash(other), Equal(other), Alloc(AllocTraits::select_on_container_copy_construction(other)),
    size_(other.size_), list_(other.list_), load_factor_(other.load_factor_),
    max_load_factor_(other.max_load_factor_), buckets_(other.buckets_.size(), list_.end()) {

    for (auto it = list_.begin(); it != list_.end(); ++it) {
      size_t hash = Hash::operator()(it->first);
      size_t bucket_index = hash % buckets_.size();
      if (buckets_[bucket_index] == list_.end()) {
        buckets_[bucket_index] = it;
      }
    }
  }

  UnorderedMap(UnorderedMap&& other) noexcept
    : Hash(std::move(other)), Equal(std::move(other)), Alloc(std::move(other)),
    size_(other.size_), list_(std::move(other.list_)), load_factor_(other.load_factor_),
    max_load_factor_(other.max_load_factor_), buckets_(std::move(other.buckets_)) {
    other.size_ = 0;
  }

  ~UnorderedMap() {}

  UnorderedMap& operator=(const UnorderedMap& other) {
    auto copy = other;
    swap(copy);
    return *this;
  }

  UnorderedMap& operator=(UnorderedMap&& other) noexcept {
    if (this == &other) {
      return *this;
    }

    Hash::operator=(std::move(other));
    Equal::operator=(std::move(other));
    if (AllocTraits::propagate_on_container_move_assignment::value) {
      Alloc::operator=(std::move(other));
    }
    size_ = other.size_;
    list_ = std::move(other.list_);
    other.size_ = 0;

    load_factor_ = other.load_factor_;
    max_load_factor_ = other.max_load_factor_;

    buckets_ = std::move(other.buckets_);
    for (size_t i = 0; i < buckets_.size(); ++i) {
      if (buckets_[i] == other.list_.end()) {
        buckets_[i] = list_.end();
      }
    }

    return *this;
  }

  Value& operator[](const Key& key) {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    iterator it = find_in_bucket(key, bucket_index);

    if (it == list_.end()) {
      return emplace(key, Value()).first->second;
    }

    return it->second;
  }

  Value& operator[](Key&& key) {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    iterator it = find_in_bucket(key, bucket_index);

    if (it == list_.end()) {
      return emplace(std::move(key), Value()).first->second;
    }

    return it->second;
  }

  Value& at(const Key& key) {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    iterator it = find_in_bucket(key, bucket_index);
    if (it == list_.end()) {
      throw std::out_of_range("Key not found");
    }
    return it->second;
  }

  const Value& at(const Key& key) const {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    const_iterator it = find_in_bucket(key, bucket_index);
    if (it == list_.end()) {
      throw std::out_of_range("Key not found");
    }
    return it->second;
  }

  size_t size() const {
    return list_.size();
  }

  iterator begin() {
    return list_.begin();
  }

  iterator end() {
    return list_.end();
  }

  const_iterator begin() const {
    return list_.begin();
  }

  const_iterator end() const {
    return list_.end();
  }

  const_iterator cbegin() const {
    return list_.cbegin();
  }

  const_iterator cend() const {
    return list_.cend();
  }

  std::pair<iterator, bool> insert(const NodeType& value) {
    size_t hash = Hash::operator()(value.first);
    size_t bucket_index = hash % buckets_.size();
    iterator it = find_in_bucket(value.first, bucket_index);
    if (it != list_.end()) {
      return { it, false };
    }

    ++size_;
    load_factor_ = static_cast<double>(size_) / buckets_.size();
    if (load_factor_ > max_load_factor_) {
      rehash(buckets_.size() * 2);
    }

    bucket_index = hash % buckets_.size();
    buckets_[bucket_index] = list_.insert(buckets_[bucket_index], value);

    return { buckets_[bucket_index], true };
  }

  std::pair<iterator, bool> insert(NodeType&& node) {
    return emplace(std::move(node));
  }

  template<typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    iterator new_node = list_.emplace(end(), std::forward<Args>(args)...);

    size_t hash = Hash::operator()(new_node->first);
    size_t bucket_index = hash % buckets_.size();
    iterator it = find_in_bucket(new_node->first, bucket_index);

    if (it != list_.end() && it != new_node) {
      list_.erase(new_node);
      return { it, false };
    }

    list_.MoveNode(new_node, buckets_[bucket_index]);
    --buckets_[bucket_index];

    ++size_;
    load_factor_ = static_cast<double>(size_) / buckets_.size();
    if (load_factor_ > max_load_factor_) {
      rehash(buckets_.size() * 2);
    }

    return { buckets_[bucket_index], true };
  }

  template<typename InputIt>
  void insert(InputIt first, InputIt last) {
    for (auto it = first; it != last; ++it) {
      insert(*it);
    }
  }

  iterator erase(const_iterator pos) {
    size_t bucket_index = Hash::operator()(pos->first) % buckets_.size();
    if (buckets_[bucket_index] == pos) {
      ++buckets_[bucket_index];

      if (buckets_[bucket_index] != end() &&
        Hash::operator()(buckets_[bucket_index]->first) % buckets_.size() != bucket_index) {
        buckets_[bucket_index] = end();
      }
    }
    --size_;
    return iterator(list_.erase(pos));
  }

  iterator erase(const_iterator first, const_iterator last) {
    for (auto it = first; it != last;) {
      it = erase(it);
    }

    return iterator(last);
  }

  iterator find(const Key& key) {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    return find_in_bucket(key, bucket_index);
  }

  const_iterator find(const Key& key) const {
    size_t hash = Hash::operator()(key);
    size_t bucket_index = hash % buckets_.size();
    return find_in_bucket(key, bucket_index);
  }

  void reserve(size_t count) {
    if (count > buckets_.size()) {
      rehash(count);
    }
  }

  double load_factor() const {
    return static_cast<double>(size_) / buckets_.size();
  }

  double max_load_factor() const {
    return max_load_factor_;
  }

  void max_load_factor(double ml) {
    max_load_factor_ = ml;
    if (load_factor() > max_load_factor_) {
      rehash(buckets_.size() * 2);
    }
  }

  void swap(UnorderedMap& other) {
    std::swap(static_cast<Hash&>(*this), static_cast<Hash&>(other));
    std::swap(static_cast<Equal&>(*this), static_cast<Equal&>(other));
    if (AllocTraits::propagate_on_container_std::swap::value) {
      std::swap(static_cast<Alloc&>(*this), static_cast<Alloc&>(other));
    }
    std::swap(size_, other.size_);
    list_.swap(other.list_);
    std::swap(load_factor_, other.load_factor_);
    std::swap(max_load_factor_, other.max_load_factor_);
    std::swap(buckets_, other.buckets_);
  }

private:
  iterator find_in_bucket(const Key& key, size_t bucket_index) const {
    const_iterator it = buckets_[bucket_index];
    while (it != list_.end() && Hash::operator()(it->first) % buckets_.size() == bucket_index) {
      if (Equal::operator()(it->first, key)) {
        return iterator(it);
      }
      ++it;
    }
    return iterator(list_.end());
  }
  

  template<typename K, typename V>
  std::pair<iterator, bool> create(K&& key, V&& value, size_t bucket_index) {
    if (load_factor() > max_load_factor_) {
      rehash(buckets_.size() * 2);
      bucket_index = Hash::operator()(key) % buckets_.size();
    }

    NodeType node(std::forward<K>(key), std::forward<V>(value));
    buckets_[bucket_index] = list_.emplace(buckets_[bucket_index], std::move(node));

    ++size_;
    return { buckets_[bucket_index], true };
  }


  UnorderedMap(size_t buckets_count, const UnorderedMap& other)
    : Hash(other), Equal(other), Alloc(other), size_(0),
    list_(other.list_.get_allocator()), load_factor_(0),
    max_load_factor_(other.max_load_factor_), buckets_(buckets_count, list_.end()) {}  //  for rehash

  void rehash(size_t new_buckets_count) {
    std::vector<std::vector<iterator>> groups(new_buckets_count);

    for (auto it = list_.begin(); it != list_.end(); ++it) {
      size_t hash = Hash::operator()(it->first);
      size_t bucket_index = hash % new_buckets_count;
      groups[bucket_index].push_back(it);
    }

    buckets_.resize(new_buckets_count);

    for (size_t i = 0; i < groups.size(); ++i) {
      if (groups[i].empty()) {
        buckets_[i] = list_.end();
        continue;
      }

      buckets_[i] = groups[i].front();
      for (auto it : groups[i]) {
        list_.MoveNode(it, list_.end());
      }
    }
  }
};
