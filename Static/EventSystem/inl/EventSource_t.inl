#ifndef ULTREALITY_UTILITIES_EVENT_SOURCE_T_INL
#define ULTREALITY_UTILITIES_EVENT_SOURCE_T_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Utilities
{
	template <CEventDispatcherCompatible EventTypeBase>
	template<typename EventType>
	FORCE_INLINE auto& EventSource_t<EventTypeBase>::GetDispatcher()
	{
		return std::get<EventDispatcher_t<EventType>>(m_eventDispatchers);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::Subscribe(void(*func)(const EventType&), Events... eventTypes)
	{
		GetDispatcher<EventType>().Subscribe(func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::Unsubscribe(void(*func)(const EventType&), Events... eventTypes)
	{
		GetDispatcher<EventType>().Unsubscribe(func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<class T, typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::Subscribe(T* instance, void(T::* func)(const EventType&), Events... eventTypes)
	{
		GetDispatcher<EventType>().Subscribe(instance, func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<class T, typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::Unsubscribe(T* instance, void(T::* func)(const EventType&), Events... eventTypes)
	{
		GetDispatcher<EventType>().Unsubscribe(instance, func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_lvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::Subscribe(Callable&& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().Subscribe(std::forward<Callable>(callable), eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_rvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::Subscribe(const Callable& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().Subscribe(callable, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_lvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::Unsubscribe(Callable&& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().Unsubscribe(std::forward<Callable>(callable), eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_rvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::Unsubscribe(const Callable& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().Unsubscribe(callable, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncSubscribe(void(*func)(const EventType& event), Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncSubscribe(func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncUnsubscribe(void(*func)(const EventType& event), Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncUnsubscribe(func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<class T, typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncSubscribe(T* instance, void(T::* func)(const EventType& event), Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncSubscribe(instance, func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<class T, typename EventType, typename... Events>
	requires (std::is_same_v<EnumType, Events> && ...)
	&& CDerivedEvent<EventTypeBase, EventType>
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncUnsubscribe(T* instance, void(T::* func)(const EventType& event), Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncUnsubscribe(instance, func, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_lvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncSubscribe(Callable&& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncSubscribe(std::forward<Callable>(callable), eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_rvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncSubscribe(const Callable& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncSubscribe(callable, eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_lvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncUnsubscribe(Callable&& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncUnsubscribe(std::forward<Callable>(callable), eventTypes...);
	}

	template <CEventDispatcherCompatible EventTypeBase>
	template<typename Callable, typename... Events>
	requires ((std::is_same_v<EnumType, Events> && ...)
	&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
	&& !std::is_rvalue_reference_v<Callable>)
	FORCE_INLINE void EventSource_t<EventTypeBase>::SyncUnsubscribe(const Callable& callable, Events... eventTypes)
	{
		GetDispatcher<EventType>().SyncUnsubscribe(callable, eventTypes...);
	}
}

#endif // !ULTREALITY_UTILITIES_EVENT_SOURCE_T_INL
