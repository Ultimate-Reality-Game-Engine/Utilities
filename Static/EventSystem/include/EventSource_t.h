#ifndef ULTREALITY_UTILITIES_EVENT_SOURCE_T_H
#define ULTREALITY_UTILITIES_EVENT_SOURCE_T_H

#include <EventDispatcher_t.h>

namespace UltReality::Utilities
{
	template <CEventDispatcherCompatible EventTypeBase>
	class EventSource_t
	{
	protected:
		/// <summary>
		/// Encapsulated event dispatcher
		/// </summary>
		EventDispatcher_t<EventTypeBase> m_eventDispatcher;

	public:
		/// <summary>
		/// Call to subscribe a free function listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void Subscribe(void(*func)(const EventType&), Events... eventTypes);

		/// <summary>
		/// Unsubscribe a free function listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void Unsubscribe(void(*func)(const EventType&), Events... eventTypes);

		/// <summary>
		/// Call to subscribe a member function listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="instance">An instance of the object with the member function to point to</param>
		/// <param name="func">A pointer to a member function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<class T, typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void Subscribe(T* instance, void(T::* func)(const EventType&), Events... eventTypes);

		/// <summary>
		/// Unsubscribe a member function listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="instance">An instance of the object with the member function pointed to</param>
		/// <param name="func">A pointer to a member function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<class T, typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void Unsubscribe(T* instance, void(T::* func)(const EventType&), Events... eventTypes);

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_lvalue_reference_v<Callable>)
		void Subscribe(Callable&& callable, Events... eventTypes);

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_rvalue_reference_v<Callable>)
		void Subscribe(const Callable& callable, Events... eventTypes);

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_lvalue_reference_v<Callable>)
		void Unsubscribe(Callable&& callable, Events... eventTypes);

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_rvalue_reference_v<Callable>)
		void Unsubscribe(const Callable& callable, Events... eventTypes);

		/// <summary>
		/// Call to subscribe a free function listener for one or many supported event types
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void SyncSubscribe(void(*func)(const EventType& event), Events... eventTypes);

		/// <summary>
		/// Unsubscribe a free function listener from specified events
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void SyncUnsubscribe(void(*func)(const EventType& event), Events... eventTypes);

		/// <summary>
		/// Call to subscribe a member function listener for one or many supported event types
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="instance">An instance of the object with the member function to point to</param>
		/// <param name="func">A pointer to a member function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<class T, typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void SyncSubscribe(T* instance, void(T::* func)(const EventType& event), Events... eventTypes);

		/// <summary>
		/// Unsubscribe a member function listener from specified events
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="instance">An instance of the object with the member function pointed to</param>
		/// <param name="func">A pointer to a member function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<class T, typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		void SyncUnsubscribe(T* instance, void(T::* func)(const EventType& event), Events... eventTypes);

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_lvalue_reference_v<Callable>)
		void SyncSubscribe(Callable&& callable, Events... eventTypes);

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_rvalue_reference_v<Callable>)
		void SyncSubscribe(const Callable& callable, Events... eventTypes);

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_lvalue_reference_v<Callable>)
		void SyncUnsubscribe(Callable&& callable, Events... eventTypes);

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires ((std::is_same_v<EnumType, Events> && ...)
		&& !std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>
		&& !std::is_rvalue_reference_v<Callable>)
		void SyncUnsubscribe(const Callable& callable, Events... eventTypes);
	};
}

#include <EventSource_t.inl>

#endif // !ULTREALITY_UTILITIES_EVENT_SOURCE_T_H
