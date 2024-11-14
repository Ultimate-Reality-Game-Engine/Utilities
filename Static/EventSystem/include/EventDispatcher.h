#ifndef ULTREALITY_UTILITIES_EVENT_DISPATCHER_H
#define ULTREALITY_UTILITIES_EVENT_DISPATCHER_H

#include <concepts>
#include <queue>
#include <unordered_map>
#include <vector>
#include <memory>
#include <mutex>
#include <algorithm>
#include <type_traits>
#include <functional>
#include <typeindex>

#include <EventBase_t.h>

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Utilities
{
	/*
		File: EventDispatcher.h

		Description:
		Defines a flexible, type-safe event dispatcher system that manages and delivers
		events to subscribed listeners in a multithreaded environment.

		How it works:
		- The 'EventDispatcher<EventTypeBase>' template class manages and dispatches events of types
		  derived from a specified 'EventBase_t' specialization, enabling a common base for various
		  event types with type-specific data.
		
		Template Parameters:
		- 'EventTypeBase': The base event type, derived from 'EventBase_t', with an associated
		  'EnumType' (enum class) representing specific event types.
		
		Constructor Parameters:
		- 'batchSize': Defines the number of events to process per dispatch batch. This number is split up
		  between the two internal queues (unsynchronized and synchronized). If one queue has less than 
		  its share of events loaded, the other queue will get that allotment for that batch dispatch.
		
		Features:
		- Queues events for either synchronous or asynchronous processing.
		- Registers listeners to handle events based on event type.
		- Supports batch processing for efficient event handling, with optional synchronization for multithreaded environments.
		
		Thread Safety:
		- Thread-safe variants of event queuing, subscribing, and unsubscribing are available through methods prefixed with 'Sync'.
		- Mutexes ensure thread-safe access to the synchronized event queue and listeners.
		- Listeners must handle concurrency independently, ensuring their callbacks are thread-safe.

		Recommendations:
		- Avoid using 'Sync' methods if the caller is on the same thread as the 'EventDispatcher' instance, as this reduces unnecessary locking overhead.
		- Use 'Sync' methods when accessing the dispatcher from different threads to prevent race conditions.

		Example Usage:
		- Define a 'PlatformEvents' enum class for OS events like 'WindowResize' and 'KeyPress'.
		- Define 'PlatformMessageEvent', which inherits from 'EventBase_t<PlatformEvents>', and use derived types like 'WindowResizeEvent' to include specific data for resize events.
		- Register listeners to handle these events by type, allowing organized and efficient event handling.

			// Example: Creating and using an EventDispatcher instance for PlatformMessageEvent type events
			EventDispatcher<PlatformMessageEvent> dispatcher;

			// Subscribe a listener for WindowResizeEvent
			dispatcher.Subscribe<WindowResizeEvent>([](const WindowResizeEvent& event) {
				// Handle window resize
			});

			// Queue a WindowResizeEvent
			dispatcher.QueueEvent(WindowResizeEvent(...));

			// Process and dispatch queued events in batches
			dispatcher.DispatchBatch();
	*/

	/// <summary>
	/// Concept to enforce that T is derived from EventBase_t type
	/// </summary>
	template<typename T>
	concept CEventBase = requires 
	{
		typename T::EnumType; // Ensures the type has an EnumType member
	} 
	&& CEnumClass<typename T::EnumType> // Ensures that the EnumType is an enum class type
	&& std::is_base_of_v<EventBase_t<typename T::EnumType>, T>; // Ensures T is derived from EventBase_t<EnumType>

	/// <summary>
	/// Concept to enforce that T is derived from EventBase_t and valid
	/// </summary>
	template<typename T>
	concept CEventDispatcherCompatible = CEventBase<T>;

	/// <summary>
	/// Concept to enforce that T is derived from the BaseType, which is compatible with EventDispatcher
	/// </summary>
	template<typename BaseType, typename T>
	concept CDerivedEvent = CEventBase<BaseType> && std::is_base_of_v<BaseType, T>;

	/// <summary>
	/// A class that handles the dispatching of event objects to subscribed listeners
	/// </summary>
	/// <typeparam name="...Events">A variadic template parameter pack of types that satisfy the <seealso cref="UltReality.Utilities.CEventType"/> concept. These are the event that manager supports</typeparam>
	template <CEventDispatcherCompatible EventTypeBase>
	class EventDispatcher
	{
	public:
		using EnumType = typename EventTypeBase::EnumType;

	private:
		struct IDelegate
		{
			virtual ~IDelegate() = default;
			virtual void operator()(const EventTypeBase&) const = 0;
			//virtual bool operator==(const IDelegate&) const = 0;
		};

		/// <summary>
		/// A wrapper around a free function pointer and utility operators to invoke as a functor
		/// </summary>
		template<typename EventType>
		requires CDerivedEvent<EventTypeBase, EventType>
		struct FreeDelegate : public IDelegate
		{
			void(*freeFunc)(const EventType&);

			/// <summary>
			/// Construct to point to a free function
			/// </summary>
			/// <param name="func">A free function pointer</param>
			explicit FreeDelegate(void(*func)(const EventType&)) : freeFunc(func) {}

			/// <summary>
			/// Will invoke the underlying function pointer with the provided event instance
			/// </summary>
			/// <param name="event">A reference to an event instance that derives from the base event type for this dispatcher</param>
			void operator()(const EventTypeBase& event) const override
			{
				freeFunc(static_cast<const EventType&>(event));
			}

			/// <summary>
			/// Compare one Delegate to another to determine if they point to the same function
			/// </summary>
			/// <param name="other">A reference to a Delegate instance to compare to</param>
			/// <returns>True if the Delegates are the same, false otherwise</returns>
			/*bool operator==(const IDelegate& other) const
			{
				return freeFunc == static_cast<const FreeDelegate&>(other).freeFunc;
			}*/
		};

		/// <summary>
		/// A wrapper around member function pointers and utility operators to invoke as a functor
		/// </summary>
		template<class T, typename EventType>
		requires CDerivedEvent<EventTypeBase, EventType>
		struct MemberDelegate : IDelegate
		{
			void(T::*memberFunc)(const EventType&);

			// Instance for member function calls
			T* m_instance;

			/// <summary>
			/// Construct to point to a member function
			/// </summary>
			/// <typeparam name="T">The type of object with the member function to point to</typeparam>
			/// <param name="instance">Instance of the specified type to call the member function on</param>
			/// <param name="func">A member function pointer</param>
			explicit MemberDelegate(T* instance, void(T::*func)(const EventType&)) : m_instance(instance), memberFunc(func) {}

			/// <summary>
			/// Will invoke the underlying function pointer with the provided event instance
			/// </summary>
			/// <param name="event">A reference to an event instance that derives from the base event type for this dispatcher</param>
			void operator()(const EventTypeBase& event) const override
			{
				(m_instance->*memberFunc)(static_cast<const EventType&>(event));
			}

			/// <summary>
			/// Compare one Delegate to another to determine if they point to the same function and reference the same object instance if applicable
			/// </summary>
			/// <param name="other">A reference to a Delegate instance to compare to</param>
			/// <returns>True if the Delegates are the same, false otherwise</returns>
			/*bool operator==(const IDelegate& other) const override
			{
				return m_instance == static_cast<const MemberDelegate&>(other).m_instance
					&& (memberFunc == static_cast<const MemberDelegate&>(other).memberFunc);
			}*/
		};

		/// <summary>
		/// A wrapper around function pointers and utility operators to invoke as a functor
		/// </summary>
		template<typename Callable, typename EventType>
		requires CDerivedEvent<EventTypeBase, EventType>
		struct Delegate_t : public IDelegate
		{
			Callable callable;

			explicit Delegate_t(Callable&& c) : callable(std::move(c)) {}

			/// <summary>
			/// Will invoke the underlying function pointer with the provided event instance
			/// </summary>
			/// <param name="event">A reference to an event instance that derives from the base event type for this dispatcher</param>
			void operator()(const EventTypeBase& event) const override
			{
				callable(static_cast<const EventType&>(event));
			}

			/// <summary>
			/// Compare one Delegate to another to determine if they point to the same function and reference the same object instance if applicable
			/// </summary>
			/// <param name="other">A reference to a Delegate instance to compare to</param>
			/// <returns>True if the Delegates are the same, false otherwise</returns>
			/*bool operator==(const IDelegate& other) const override
			{
				if (const auto* otherDelegate = dynamic_cast<const Delegate_t*>(&other))
				{
					return &callable == &otherDelegate->callable;
				}
				return false;
			}*/
		};

		template<typename T>
		struct function_traits;

		// Specialization fot std::function
		template<typename Ret, typename... Args>
		struct function_traits<std::function<Ret(Args...)>>
		{
			using result_type = Ret;
			using args_tuple = std::tuple<Args...>;
		};

		// Fallback for callable objects, including lambdas
		template<typename Callable>
		struct function_traits : function_traits<decltype(&Callable::operator())> {};

		// Specialization for lambdas or any callable objects with 'operator()'
		template<typename Ret, typename ClassType, typename... Args>
		struct function_traits<Ret(ClassType::*)(Args...) const>
		{
			using result_type = Ret;
			using args_tuple = std::tuple<Args...>;
		};

		// Helper to get the type of the first argument (event type in the case of lambdas passed to Subscribe)
		template<typename Callable>
		using first_argument_t = std::tuple_element_t<0, typename function_traits<Callable>::args_tuple>;


		size_t m_batchSize; // The maximum number of events to process in one batch

		// Unsynchronized collections for same-thread access
		std::queue<std::shared_ptr<EventTypeBase>> m_eventQueue; // A queue of event instances that this dispatcher will process
		std::unordered_map<EnumType, std::unordered_map<size_t, std::shared_ptr<IDelegate>>> m_listeners; // A map of listeners associated with the event types that will trigger them

		// Synchronized collections for cross-thread access
		std::queue<std::shared_ptr<EventTypeBase>> m_syncEventQueue;
		std::unordered_map<EnumType, std::unordered_map<size_t, std::shared_ptr<IDelegate>>> m_syncListeners;

		// Mutexes for synchronized collections
		mutable std::mutex m_eventQueueMutex;
		mutable std::mutex m_listenerMutex;

	public:
		EventDispatcher(size_t batchSize = 60) : m_batchSize(batchSize)
		{}

		/// <summary>
		/// Set the batch size. It is recommended to set a multiple of 2, as there are two internal queues that each get allocated half this number ber batch dispatch
		/// </summary>
		/// <param name="batchSize">The number of events to dispatch per batch in total</param>
		void SetBatchSize(size_t batchSize)
		{
			m_batchSize = batchSize;
		}

		/// <summary>
		/// Check to see if total queue volume is 0
		/// </summary>
		/// <returns>True if there are no queued events, false otherwise</returns>
		bool IsQueueEmpty()
		{
			std::lock_guard<std::mutex> lock(m_eventQueueMutex);
			return m_eventQueue.empty() && m_syncEventQueue.empty();
		}

		/// <summary>
		/// Get the total queue volume
		/// </summary>
		/// <returns>The total number of events queued for dispatch</returns>
		size_t QueueSize()
		{
			std::lock_guard<std::mutex> lock(m_eventQueueMutex);
			return m_eventQueue.size() + m_syncEventQueue.size();
		}

		/// <summary>
		/// Call to add an event to the queue for the manager to dispatch to subscribers
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="event">An instance of a supported event type</param>
		template<typename T>
		requires CDerivedEvent<EventTypeBase, T>
		FORCE_INLINE void QueueEvent(const T& event)
		{
			m_eventQueue.emplace(std::make_shared<T>(event));
		}

		/// <summary>
		/// Call to add an event to the queue for the manager to dispatch to subscribers
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="event">An instance of a supported event type</param>
		template<typename T>
		requires CDerivedEvent<EventTypeBase, T>
		FORCE_INLINE void SyncQueueEvent(const T& event)
		{
			std::lock_guard<std::mutex> lock(m_eventQueueMutex);
			m_syncEventQueue.emplace(std::make_shared<T>(event));
		}

		/// <summary>
		/// Call to subscribe a free function listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		FORCE_INLINE void Subscribe(void(*func)(const EventType&), Events... eventTypes)
		{
			((m_listeners[eventTypes][Hash(func)] = std::make_shared<FreeDelegate<EventType>>(func)), ...);
		}

		/// <summary>
		/// Unsubscribe a free function listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		FORCE_INLINE void Unsubscribe(void(*func)(const EventType&), Events... eventTypes)
		{
			(m_listeners[eventTypes].erase(Hash(func)), ...);
		}

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
		FORCE_INLINE void Subscribe(T* instance, void(T::*func)(const EventType&), Events... eventTypes)
		{
			((m_listeners[eventTypes][Hash(instance, func)] = std::make_shared<MemberDelegate<T, EventType>>(instance, func)), ...);
		}

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
		FORCE_INLINE void Unsubscribe(T* instance, void(T::*func)(const EventType&), Events... eventTypes)
		{
			(m_listeners[eventTypes].erase(Hash(instance, func)), ...);
		}

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		Subscribe(Callable&& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<Callable>>;

			((m_listeners[eventTypes][Hash(callable)] = std::make_shared<Delegate_t<Callable, EventType>>(std::forward<Callable>(callable))), ...);
		}

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		Subscribe(const Callable& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<std::remove_reference_t<decltype(Callable)>>>;

			((m_listeners[eventTypes][Hash(callable)] = std::make_shared<Delegate_t<std::remove_reference_t<decltype(Callable)>, EventType>>(callable)), ...);
		}

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		Unsubscribe(Callable&& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<Callable>>;

			(m_listeners[eventTypes].erase(Hash(callable)), ...);
		}

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		Unsubscribe(const Callable& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<std::decay_t<decltype(Callable)>>>;

			(m_listeners[eventTypes].erase(Hash(callable)), ...);
		}

		/// <summary>
		/// Call to subscribe a free function listener for one or many supported event types
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		FORCE_INLINE void SyncSubscribe(void(*func)(const EventType& event), Events... eventTypes)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			((m_syncListeners[eventTypes][Hash(func)] = std::make_shared<FreeDelegate<EventType>>(func)), ...);
		}

		/// <summary>
		/// Unsubscribe a free function listener from specified events
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="func">A pointer to a free function that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename EventType, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		&& CDerivedEvent<EventTypeBase, EventType>
		FORCE_INLINE void SyncUnsubscribe(void(*func)(const EventType& event), Events... eventTypes)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(m_syncListeners[eventTypes].erase(Hash(func)), ...);
		}

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
		FORCE_INLINE void SyncSubscribe(T* instance, void(T::*func)(const EventType& event), Events... eventTypes)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);

			((m_syncListeners[eventTypes][Hash(instance, func)] = std::make_shared<MemberDelegate<T, EventType>>(instance, func)), ...);
		}

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
		FORCE_INLINE void SyncUnsubscribe(T* instance, void(T::*func)(const EventType& event), Events... eventTypes)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);

			(m_syncListeners[eventTypes].erase(Hash(instance, func)), ...);
		}

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		SyncSubscribe(Callable&& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<Callable>>;

			std::lock_guard<std::mutex> lock(m_listenerMutex);
			((m_syncListeners[eventTypes][Hash(callable)] = std::make_shared<Delegate_t<Callable, EventType>>(std::forward<Callable>(callable))), ...);
		}

		/// <summary>
		/// Call to subscribe a general callable object (lambda, std::function) listener for one or many supported event types
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to subscribe to</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		SyncSubscribe(const Callable& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<std::decay_t<decltype(Callable)>>>;

			std::lock_guard<std::mutex> lock(m_listenerMutex);
			((m_syncListeners[eventTypes][Hash(callable)] = std::make_shared<Delegate_t<std::decay_t<decltype(Callable)>, EventType>>(callable)), ...);
		}

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		SyncUnsubscribe(Callable&& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<Callable>>;

			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(m_syncListeners[eventTypes].erase(Hash(callable)), ...);
		}

		/// <summary>
		/// Unsubscribe a general callable object (lambda, std::function) listener from specified events
		/// Thread safe. Can be called by concurrently operating threads
		/// </summary>
		/// <param name="callable">An rvalue reference to a callable object that takes an 'EventTypeBase' compatible const reference</param>
		/// <param name="eventTypes">A pack of 'EnumType' values specifying the event types to unsubscribe from</param>
		template<typename Callable, typename... Events>
		requires (std::is_same_v<EnumType, Events> && ...)
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, void>
		SyncUnsubscribe(const Callable& callable, Events... eventTypes)
		{
			using EventType = std::remove_cvref_t<first_argument_t<std::decay_t<decltype(Callable)>>>;

			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(m_syncListeners[eventTypes].erase(Hash(callable)), ...);
		}

		/// <summary>
		/// Dispatches a batch of events to subscribed listeners
		/// <seealso cref="batchSize"/> events or <seealso cref="m_eventQueue"/>.size(), whichever is less
		/// </summary>
		void DispatchBatch()
		{
			std::lock_guard<std::mutex> queueLock(m_eventQueueMutex);
			std::lock_guard<std::mutex> listenerLock(m_listenerMutex);

			size_t halfBatchSize = m_batchSize / 2;

			size_t count = std::min(m_eventQueue.size(), halfBatchSize);
			size_t syncCount = std::min(m_syncEventQueue.size(), halfBatchSize);

			size_t remaining = m_batchSize - (count + syncCount);

			// If there's still room in the batch, take the remainder from the other queue
			if (remaining > 0)
			{
				if (count < halfBatchSize)
				{
					syncCount = std::min(syncCount + remaining, m_syncEventQueue.size());
				}
				else
				{
					count = std::min(count + remaining, m_eventQueue.size());
				}
			}

			for (size_t i = 0; i < count; i++)
			{
				const auto& event = m_eventQueue.front();
				NotifySubscribers(event, m_listeners);
				NotifySubscribers(event, m_syncListeners);
				m_eventQueue.pop();
			}

			for (size_t i = 0; i < syncCount; i++)
			{
				const auto& event = m_syncEventQueue.front();
				NotifySubscribers(event, m_listeners);
				NotifySubscribers(event, m_syncListeners);
				m_syncEventQueue.pop();
			}
		}

		/// <summary>
		/// Dispatch all events that are queued to subscribed listeners
		/// </summary>
		void DispatchAll()
		{
			std::lock_guard<std::mutex> queueLock(m_eventQueueMutex);
			std::lock_guard<std::mutex> listenerLock(m_listenerMutex);

			while (!m_eventQueue.empty())
			{
				const auto& event = m_eventQueue.front();
				NotifySubscribers(event, m_listeners);
				m_eventQueue.pop();
			}

			while (!m_syncEventQueue.empty())
			{
				const auto& event = m_syncEventQueue.front();
				NotifySubscribers(event, m_syncListeners);
				m_syncEventQueue.pop();
			}
		}

	private:
		/// <summary>
		/// Helper method to call the callback method's for the subscribed listeners associated with the event
		/// </summary>
		/// <param name="event">An event instance being processed</param>
		FORCE_INLINE void NotifySubscribers(const std::shared_ptr<EventTypeBase>& event, 
			const std::unordered_map<EnumType, std::unordered_map<size_t, std::shared_ptr<IDelegate>>>& listeners)
		{
			auto it = listeners.find(event->type);
			if (it != listeners.end())
			{
				for (const auto& callback : it->second)
				{
					callback.second->operator()(*event);
				}
			}
		}

		/// <summary>
		/// Helper method for removing a Delegate from one of the vectors in the maps
		/// </summary>
		/// <param name="listeners">Reference to the vector to search and remove from</param>
		/// <param name="delegate">Instance to remove, used to compare to instance in vector</param>
		/*FORCE_INLINE void RemoveListener(std::unordered_map<size_t, std::shared_ptr<IDelegate>>& listeners, size_t hash)
		{
			listeners.erase(hash);
		}*/

		template<typename EventType>
		FORCE_INLINE size_t Hash(void(*func)(const EventType&))
		{
			return std::hash<void*>()(reinterpret_cast<void*>(func));
		}

		template<typename T, typename EventType>
		FORCE_INLINE size_t Hash(T* instance, void(T::*func)(const EventType&))
		{
			size_t instanceHash = std::hash<T*>()(instance);
			size_t funcHash = std::type_index(typeid(func)).hash_code();
			return instanceHash ^ (funcHash << 1);
		}

		template<typename Callable>
		FORCE_INLINE typename std::enable_if_t<!std::is_function_v<std::remove_pointer_t<std::decay_t<Callable>>>, size_t>
		Hash(const Callable& callable)
		{
			return typeid(callable).hash_code();
		}
	};
}
#endif // !ULTREALITY_UTILITIES_EVENT_DISPATCHER_H
