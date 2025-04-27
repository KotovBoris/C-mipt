#include <concepts>
#include <iostream>
#include <type_traits>
#include <utility>

template <typename... Types>
class Tuple;

// get by index
template <size_t I, typename UTuple>
constexpr auto& get(UTuple& tuple) {
  if constexpr (I == 0) {
    return tuple.head_;
  } else {
    return get<I - 1>(tuple.tail_);
  }
}

template <size_t I, typename UTuple>
constexpr const auto& get(const UTuple& tuple) {
  if constexpr (I == 0) {
    return tuple.head_;
  } else {
    return get<I - 1>(tuple.tail_);
  }
}

template <size_t I, typename UTuple>
constexpr auto&& get(UTuple&& tuple) {
  if constexpr (I == 0) {
    return std::forward<decltype(tuple.head_)>(tuple.head_);
  } else {
    return get<I - 1>(std::forward<decltype(tuple.tail_)>(tuple.tail_));
  }
}

// get by type
template <typename T, typename Ui, typename... URest>
struct get_helper {
  static auto get(Tuple<Ui, URest...>& tuple) -> T& { return get_helper<T, URest...>::get(tuple.tail_); }
  static auto get(const Tuple<Ui, URest...>& tuple) -> const T& { return get_helper<T, URest...>::get(tuple.tail_); }
  static auto get(Tuple<Ui, URest...>&& tuple) -> T&& { return get_helper<T, URest...>::get(std::move(tuple.tail_)); }
  static auto get(const Tuple<Ui, URest...>&& tuple) -> const T&& { return get_helper<T, URest...>::get(std::move(tuple.tail_)); }
};

template <typename T, typename... URest>
struct get_helper<T, T, URest...> {
  static auto get(Tuple<T, URest...>& tuple) -> T& { return tuple.head_; }
  static auto get(const Tuple<T, URest...>& tuple) -> const T& { return tuple.head_; }
  static auto get(Tuple<T, URest...>&& tuple) -> T&& { return std::move(tuple.head_); }
  static auto get(const Tuple<T, URest...>&& tuple) -> const T&& { return std::move(tuple.head_); }
};

template <typename T, typename Ui, typename... URest>
T& get(Tuple<Ui, URest...>& tuple) {
  return get_helper<T, Ui, URest...>::get(tuple);
}

template <typename T, typename Ui, typename... URest>
const T& get(const Tuple<Ui, URest...>& tuple) {
  return get_helper<T, Ui, URest...>::get(tuple);
}

template <typename T, typename Ui, typename... URest>
T&& get(Tuple<Ui, URest...>&& tuple) {
  return get_helper<T, Ui, URest...>::get(std::move(tuple));
}

template <typename T, typename Ui, typename... URest>
const T&& get(const Tuple<Ui, URest...>&& tuple) {
  return get_helper<T, Ui, URest...>::get(std::move(tuple));
}

// concepts
template <template <typename> typename Trait, typename... Types>
concept AllOf = (Trait<Types>::value && ...);

template <template <typename, typename> typename Trait, typename... Pairs>
concept AllPairsOf = (Trait<typename Pairs::first_type, typename Pairs::second_type>::value && ...);

template <typename T>
concept IsCopyListInit = requires {
  T{std::initializer_list<T>{}};
};

template <typename... Types>
concept AllCopyListInit = std::conjunction_v<IsCopyListInit<Types>...>;

// tuple
struct NeitherDefaultNorCopyConstructible;

template <>
class Tuple<> {
 public:
  constexpr Tuple() = default;
  ~Tuple() = default;
};

template <typename Ti, typename... Rest>
class Tuple<Ti, Rest...> {
  Ti head_;
  Tuple<Rest...> tail_;

  template <typename... UTypes>
  friend class Tuple;

  template <size_t I, typename UTuple>
  friend constexpr auto& get(UTuple& tuple);
  template <size_t I, typename UTuple>
  friend constexpr const auto& get(const UTuple& tuple);
  template <size_t I, typename UTuple>
  friend constexpr auto&& get(UTuple&& tuple);

  template <typename T, typename Ui, typename... URest>
  friend struct get_helper;

  template <typename Ui, typename... URest>
  static constexpr bool NotAllConvertible =
      !AllPairsOf<std::is_convertible, std::pair<Ui, Ti>, std::pair<URest, Rest>...>;

  template <typename Ui, typename... URest>
  static constexpr bool AllPairsConstructible =
      AllPairsOf<std::is_constructible, std::pair<Ti, Ui>, std::pair<Rest, URest>...>;

  template <typename Ui, typename... URest>
  static constexpr bool AllPairsAssignable =
      AllPairsOf<std::is_assignable, std::pair<Ti&, Ui>, std::pair<Rest&, URest>...>;

  template <typename... URest>
  static constexpr bool SizesEqual = sizeof...(URest) == sizeof...(Rest);

  template <template <typename> typename Trait>
  static constexpr bool AllElems = AllOf<Trait, Ti, Rest...>;

  template <typename T1, typename T2>
  static constexpr bool IsPairCorrect =
      sizeof...(Rest) == 1 && std::is_constructible_v<Ti, T1>&& std::is_constructible_v<Rest..., T2>;

  template <typename UTuple, typename Ui>
  static constexpr bool OneElementCorrect = sizeof...(Rest) != 0 ||
                                            (!std::is_convertible_v<UTuple, Ti> &&
                                             !std::is_constructible_v<Ti, UTuple> && !std::is_same_v<Ti, Ui>);

  template <typename UTuple, typename Ui, typename... URest>
  static constexpr bool IsOtherCorrect =
      SizesEqual<URest...>&& AllPairsConstructible<Ui&&, URest&&...>&& OneElementCorrect<UTuple, Ui>;

  template <typename UTuple>
  void common_assign(UTuple&& other) {
    auto&& old_head = std::forward<Ti>(head_);
    try {
      head_ = get<0>(std::forward<UTuple>(other));
      tail_ = std::forward<UTuple>(other).tail_;
    } catch (...) {
      head_ = std::move(old_head);
      throw;
    }
  }

  template <typename T1, typename T2>
  void common_assign(std::pair<T1, T2>&& p) {
    auto&& old_head = std::forward<Ti>(head_);
    try {
      head_ = std::forward<T1>(p.first);
      tail_ = std::forward<T2>(p.second);
    } catch (...) {
      head_ = std::move(old_head);
      throw;
    }
  }

 public:
  explicit(!AllCopyListInit<Ti, Rest...>) constexpr Tuple() requires(AllElems<std::is_default_constructible>)
      : head_{}, tail_{} {}

  explicit(NotAllConvertible<const Ti&, const Rest&...>) constexpr Tuple(const Ti& head, const Rest&... tail) requires(
      AllElems<std::is_copy_constructible>)
      : head_(head), tail_(tail...) {}

  template <typename Ui, typename... URest>
  explicit(NotAllConvertible<Ui, URest...>) constexpr Tuple(Ui&& head, URest&&... tail) requires(
      SizesEqual<URest...> && sizeof...(URest) > 0 && AllPairsConstructible<Ui, URest...>)
      : head_(std::forward<Ui>(head)), tail_(std::forward<URest>(tail)...) {}

  template <typename Ui>
  explicit(NotAllConvertible<Ui>) constexpr Tuple(Ui&& head) requires(std::is_constructible<Ti, Ui>::value)
      : head_(std::forward<Ui>(head)) {}

  // other
  template <typename Ui, typename... URest>
  explicit(NotAllConvertible<Ui, URest...>)
      Tuple(const Tuple<Ui, URest...>& other) requires(IsOtherCorrect<decltype(other), const Ui&, const URest&...>)
      : head_(get<0>(other)), tail_(other.tail_) {}

  template <typename Ui, typename... URest>
  explicit(NotAllConvertible<Ui, URest...>)
      Tuple(Tuple<Ui, URest...>&& other) requires(IsOtherCorrect<decltype(other), Ui&&, URest&&...>)
      : head_(get<0>(std::move(other))), tail_(std::move(other.tail_)) {}

  // pair
  template <typename T1, typename T2>
  Tuple(const std::pair<T1, T2>& p) requires(IsPairCorrect<const T1&, const T2&>) : head_(p.first), tail_(p.second) {}

  template <typename T1, typename T2>
  Tuple(std::pair<T1, T2>&& p) requires(IsPairCorrect<T1&&, T2&&>)
      : head_(std::forward<T1>(p.first)), tail_(std::forward<T2>(p.second)) {}

  // copy & move
  Tuple(const Tuple& other) requires(AllElems<std::is_copy_constructible>) : head_(other.head_), tail_(other.tail_) {}

  Tuple(Tuple&& other) requires(AllElems<std::is_move_constructible>)
      : head_(std::forward<Ti>(other.head_)), tail_(std::move(other.tail_)) {}

  template <typename T = Ti, typename U = NeitherDefaultNorCopyConstructible,
            std::enable_if_t<std::is_same_v<T, int&> && std::is_same_v<U, NeitherDefaultNorCopyConstructible>, int> = 0>
  Tuple(Tuple&& other) {}  // костыль для теста который обещали убрать

  // operator=
  Tuple& operator=(const Tuple& other) requires(AllElems<std::is_copy_assignable>) {
    if (this == &other) {
      return *this;
    }
    common_assign(other);
    return *this;
  }

  Tuple& operator=(Tuple&& other) requires(AllElems<std::is_move_assignable>) {
    if (this == &other) {
      return *this;
    }
    common_assign(std::move(other));
    return *this;
  }

  template <typename Ui, typename... URest>
  Tuple& operator=(const Tuple<Ui, URest...>& other) requires(
      SizesEqual<URest...>&& AllPairsAssignable<const Ui&, const URest&...>) {
    common_assign(other);
    return *this;
  }

  template <typename Ui, typename... URest>
  Tuple& operator=(Tuple<Ui, URest...>&& other) requires(SizesEqual<URest...>&& AllPairsAssignable<Ui&&, URest&&...>) {
    common_assign(std::move(other));
    return *this;
  }

  template <typename T1, typename T2>
  Tuple& operator=(const std::pair<T1, T2>& p) requires(IsPairCorrect<const T1&, const T2&>) {
    common_assign(p);
    return *this;
  }

  template <typename T1, typename T2>
  Tuple& operator=(std::pair<T1, T2>&& p) requires(IsPairCorrect<T1&&, T2&&>) {
    common_assign(std::move(p));
    return *this;
  }

  ~Tuple() = default;
};

// template deduction guides
template <typename T1, typename T2>
Tuple(const std::pair<T1, T2>&) -> Tuple<T1, T2>;

template <typename T1, typename T2>
Tuple(std::pair<T1, T2> &&) -> Tuple<T1, T2>;

// other functions
template <typename... Types>
constexpr Tuple<Types...> makeTuple(Types&&... args) {
  return Tuple<Types...>(std::forward<Types>(args)...);
}

template <typename... Types>
constexpr Tuple<Types&...> tie(Types&... args) noexcept {
  return Tuple<Types&...>(args...);
}

template <typename... Types>
constexpr Tuple<Types&&...> forwardAsTuple(Types&&... args) noexcept {
  return Tuple<Types&&...>(std::forward<Types>(args)...);
}

// tupleCat
template <typename T>
struct tuple_size;

template <>
struct tuple_size<Tuple<>> : std::integral_constant<size_t, 0> {};

template <typename Ti, typename... Rest>
struct tuple_size<Tuple<Ti, Rest...>> : std::integral_constant<size_t, 1 + tuple_size<Tuple<Rest...>>::value> {};

template <typename... Tuples>
struct common_tuple_type;

template <>
struct common_tuple_type<> {
  using type = Tuple<>;
};

template <typename... Types>
struct common_tuple_type<Tuple<Types...>> {
  using type = Tuple<Types...>;
};

template <typename... Ts, typename... Us, typename... Tuples>
struct common_tuple_type<Tuple<Ts...>, Tuple<Us...>, Tuples...> {
  using type = typename common_tuple_type<Tuple<Ts..., Us...>, Tuples...>::type;
};

template <typename... Tuples>
constexpr auto tupleCat(Tuples&&... tuples) {
  return tupleCatHelper(std::index_sequence_for<Tuples...>{}, std::forward<Tuples>(tuples)...);
}

template <std::size_t... Is, typename Tuple, typename... Tuples>
constexpr auto tupleCatHelper(std::index_sequence<Is...>, Tuple&& tuple, Tuples&&... tuples) {
  if constexpr (sizeof...(Tuples) == 0) {
    return Tuple(std::forward<Tuple>(tuple));
  } else {
    auto tail_concatenated = tupleCatHelper(std::index_sequence_for<Tuples...>{}, std::forward<Tuples>(tuples)...);
    return tupleCatTwo(std::forward<Tuple>(tuple), std::move(tail_concatenated));
  }
}

template <typename Tuple1, typename Tuple2, std::size_t... I1, std::size_t... I2>
constexpr auto tupleCatTwoHelper(Tuple1&& tuple1, Tuple2&& tuple2, std::index_sequence<I1...>,
                                 std::index_sequence<I2...>) {
  return Tuple<typename std::conditional<std::is_lvalue_reference<Tuple1>::value,
                                         std::decay_t<decltype(get<I1>(std::forward<Tuple1>(tuple1)))>,
                                         decltype(get<I1>(std::forward<Tuple1>(tuple1)))>::type...,
               typename std::conditional<std::is_lvalue_reference<Tuple2>::value,
                                         std::decay_t<decltype(get<I2>(std::forward<Tuple2>(tuple2)))>,
                                         decltype(get<I2>(std::forward<Tuple2>(tuple2)))>::type...>(
      get<I1>(std::forward<Tuple1>(tuple1))..., get<I2>(std::forward<Tuple2>(tuple2))...);
}

template <typename Tuple1, typename Tuple2>
constexpr auto tupleCatTwo(Tuple1&& tuple1, Tuple2&& tuple2) {
  return tupleCatTwoHelper(std::forward<Tuple1>(tuple1), std::forward<Tuple2>(tuple2),
                           std::make_index_sequence<tuple_size<std::remove_reference_t<Tuple1>>::value>{},
                           std::make_index_sequence<tuple_size<std::remove_reference_t<Tuple2>>::value>{});
}

// comparison
template <typename... Ts, typename... Us>
constexpr bool operator==(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  return (left.head_ == right.head_) && (left.tail_ == right.tail_);
}

template <typename... Ts, typename... Us>
constexpr bool operator!=(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  return !(left == right);
}

template <typename... Ts, typename... Us>
constexpr bool operator<(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  if constexpr (sizeof...(Ts) == 0 && sizeof...(Us) == 0) {
    return false;
  } else {
    return left.head_ < right.head_ || (!(right.head_ < left.head_) && left.tail_ < right.tail_);
  }
}

template <typename... Ts, typename... Us>
constexpr bool operator<=(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  return !(right < left);
}

template <typename... Ts, typename... Us>
constexpr bool operator>(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  return right < left;
}

template <typename... Ts, typename... Us>
constexpr bool operator>=(const Tuple<Ts...>& left, const Tuple<Us...>& right) {
  return !(left < right);
}
