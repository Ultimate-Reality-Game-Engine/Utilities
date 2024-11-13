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
		- 'batchSize': Defines the number of events to process per dispatch batch from both
		   synchronized and unsynchronized maps (in effect the total is tice this number);
		   defaults to 60 (120 total per batch at maximum)
		
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
	/// A class that handles the dispatching of event objects to subscribed listeners
	/// </summary>
	/// <typeparam name="...Events">A variadic template parameter pack of types that satisfy the <seealso cref="UltReality.Utilities.CEventType"/> concept. These are the event that manager supports</typeparam>
	template <CEventDispatcherCompatible EventTypeBase>
	class EventDispatcher
	{
	public:
		using FreeFunction = void(*)(const EventTypeBase&);

		template<typename T>
		using MemberFunction = void(T::*)(const EventTypeBase&);

		using EnumType = typename EventTypeBase::EnumType;

	private:
		/// <summary>
		/// A wrapper around function pointers and utility operators to invoke as a functor
		/// </summary>
		struct Delegate
		{
			// Union that holds either a free function or a member function
			union
			{
				void(*freeFunc)(const EventTypeBase&);
				void(*memberFunc)(void*, const EventTypeBase&);
			};

			// Instance for member function calls
			void* m_instance = nullptr;

			/// <summary>
			/// Construct to point to a free function
			/// </summary>
			/// <param name="func">A free function pointer</param>
			Delegate(void(*func)(const EventTypeBase&)) : freeFunc(func) {}

			/// <summary>
			/// Construct to point to a member function
			/// </summary>
			/// <typeparam name="T">The type of object with the member function to point to</typeparam>
			/// <param name="instance">Instance of the specified type to call the member function on</param>
			/// <param name="func">A member function pointer</param>
			template<class T>
			Delegate(T* instance, void(T::* func)(const EventTypeBase&)) : m_instance(instance)
			{
				memberFunc = [](void* obj, const EventTypeBase& event) {
					(static_cast<T*>(obj)->*func)(event);
				};
			}

			/// <summary>
			/// Will invoke the underlying function pointer with the provided event instance
			/// </summary>
			/// <param name="event">A reference to an event instance that derives from the base event type for this dispatcher</param>
			void operator()(const EventTypeBase& event) const
			{
				if (!m_instance)
				{
					freeFunc(event);
				}
				else
				{
					memberFunc(m_instance, event);
				}
			}

			/// <summary>
			/// Compare one Delegate to another to determine if they point to the same function and reference the same object instance if applicable
			/// </summary>
			/// <param name="other">A reference to a Delegate instance to compare to</param>
			/// <returns>True if the Delegates are the same, false otherwise</returns>
			bool operator==(const Delegate& other) const
			{
				return m_instance == other.m_instance
					&& ((!m_instance && freeFunc == other.freeFunc) || (m_instance && memberFunc == other.memberFunc));
			}
		};


		size_t m_batchSize; // The maximum number of events to process in one batch

		// Unsynchronized collections for same-thread access
		std::queue<std::shared_ptr<EventTypeBase>> m_eventQueue; // A queue of event instances that this dispatcher will process
		std::unordered_map<EnumType, std::vector<Delegate>> m_listeners; // A map of listeners associated with the event types that will trigger them

		// Synchronized collections for cross-thread access
		std::queue<std::shared_ptr<EventTypeBase>> m_syncEventQueue;
		std::unordered_map<EnumType, std::vector<Delegate>> m_syncListeners;

		// Mutexes for synchronized collections
		mutable std::mutex m_eventQueueMutex;
		mutable std::mutex m_listenerMutex;

	public:
		EventDispatcher(size_t batchSize = 60) : m_batchSize(batchSize)
		{}

		/// <summary>
		/// Call to add an event to the queue for the manager to dispatch to subscribers
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="event">An instance of a supported event type</param>
		FORCE_INLINE void QueueEvent(const EventTypeBase& event)
		{
			m_eventQueue.emplace(std::make_shared<EventTypeBase>(event));
		}

		/// <summary>
		/// Call to add an event to the queue for the manager to dispatch to subscribers
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="event">An instance of a supported event type</param>
		FORCE_INLINE void SyncQueueEvent(const EventTypeBase& event)
		{
			std::lock_guard<std::mutex> lock(m_eventQueueMutex);
			m_syncEventQueue.emplace(std::make_shared<EventTypeBase>(event));
		}

		/// <summary>
		/// Call to subscribe a listener callback method for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="listener">An instance of a delegate that takes an 'EventTypeBase' compatible const reference</param>
		template<EnumType... EventTypes>
		FORCE_INLINE void Subscribe(FreeFunction listener)
		{
			(m_listeners[EventTypes].emplace_back(listener), ...);
		}

		/// <summary>
		/// Unsubscribe a listener from specified events
		/// </summary>
		template<EnumType... EventTypes>
		FORCE_INLINE void Unsubscribe(FreeFunction listener)
		{
			(RemoveListener(m_listeners[EventTypes], Delegate(listener)), ...);
		}

		/// <summary>
		/// Call to subscribe a listener callback method for one or many supported event types
		/// Not thread safe. Only call from same thread in series
		/// </summary>
		/// <param name="instance">An instance of the object with the member function to point to</param>
		/// <param name="listener">An instance of a delegate that takes an 'EventTypeBase' compatible const reference</param>
		template<class T, EnumType... EventTypes>
		FORCE_INLINE void Subscribe(T* instance, MemberFunction<T> listener)
		{
			(m_listeners[EventTypes].emplace_back(instance, listener), ...);
		}

		/// <summary>
		/// Unsubscribe a listener from specified events
		/// </summary>
		template<class T, EnumType... EventTypes>
		FORCE_INLINE void Unsubscribe(T* instance, MemberFunction<T> listener)
		{
			(RemoveListener(m_listeners[EventTypes], Delegate(instance, listener)), ...);
		}

		template<EnumType EventType, typename Callable>
		FORCE_INLINE void Subscribe(const Callable& lambda)
		{
			m_listeners[EventType].emplace_back(&lambda, &Callable::operator());
		}

		/// <summary>
		/// Call to subscribe a listener callback method for one or many supported event types
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="listener">An instance of a delegate that takes an 'EventTypeBase' compatible const reference</param>
		template<EnumType... EventTypes>
		FORCE_INLINE void SyncSubscribe(FreeFunction listener)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(m_syncListeners[EventTypes].emplace_back(listener), ...);
		}

		/// <summary>
		/// Unsubscribe a listener from specified events in thread safe manner
		/// </summary>
		template<EnumType... EventTypes>
		FORCE_INLINE void SyncUnsubscribe(FreeFunction listener)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(RemoveListener(m_syncListeners[EventTypes], Delegate(listener)), ...);
		}

		/// <summary>
		/// Call to subscribe a listener callback method for one or many supported event types
		/// Thread safe. Can be called from concurrent threads
		/// </summary>
		/// <param name="instance">An instance of the object with the member function to point to</param>
		/// <param name="listener">An instance of a delegate that takes an 'EventTypeBase' compatible const reference</param>
		template<class T, EnumType... EventTypes>
		FORCE_INLINE void SyncSubscribe(T* instance, MemberFunction<T> listener)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(m_syncListeners[EventTypes].emplace_back(instance, listener), ...);
		}

		/// <summary>
		/// Unsubscribe a listener from specified events in a thread safe manner
		/// </summary>
		template<class T, EnumType... EventTypes>
		FORCE_INLINE void SyncUnsubscribe(T* instance, MemberFunction<T> listener)
		{
			std::lock_guard<std::mutex> lock(m_listenerMutex);
			(RemoveListener(m_syncListeners[EventTypes], Delegate(instance, listener)), ...);
		}

		/// <summary>
		/// Dispatches a batch of events to subscribed listeners
		/// <seealso cref="batchSize"/> events or <seealso cref="m_eventQueue"/>.size(), whichever is less
		/// </summary>
		void DispatchBatch()
		{
			size_t count = std::min(m_eventQueue.size(), m_batchSize);

			for (size_t i = 0; i < count; i++)
			{
				const auto& event = m_eventQueue.front();
				m_eventQueue.pop();
				NotifySubscribers(event, m_listeners);
			}

			std::lock_guard<std::mutex> queueLock(m_eventQueueMutex);
			size_t syncCount = std::min(m_syncEventQueue.size(), m_batchSize);

			std::lock_guard<std::mutex> listenerLock(m_listenerMutex);
			for (size_t i = 0; i < count; i++)
			{
				const auto& event = m_syncEventQueue.front();
				m_syncEventQueue.pop();
				NotifySubscribers(event, m_syncListeners);
			}
		}

		/// <summary>
		/// Dispatch all events that are queued to subscribed listeners
		/// </summary>
		void DispatchAll()
		{
			while (!m_eventQueue.empty())
			{
				const auto& event = m_eventQueue.front();
				m_eventQueue.pop();
				NotifySubscribers(event, m_listeners);
			}

			std::lock_guard<std::mutex> queueLock(m_eventQueueMutex);
			std::lock_guard<std::mutex> listenerLock(m_listenerMutex);
			while (!m_syncEventQueue.empty())
			{
				const auto& event = m_syncEventQueue.front();
				m_syncEventQueue.pop();
				NotifySubscribers(event, m_syncListeners);
			}
		}

	private:
		/// <summary>
		/// Helper method to call the callback method's for the subscribed listeners associated with the event
		/// </summary>
		/// <param name="event">An event instance being processed</param>
		FORCE_INLINE void NotifySubscribers(const std::shared_ptr<EventTypeBase>& event, 
			const std::unordered_map<typename EventTypeBase::EnumType, std::vector<Delegate>>& listeners)
		{
			auto it = listeners.find(event->type);
			if (it != listeners.end())
			{
				for (const auto& callback : it->second)
				{
					callback(*event);
				}
			}
		}

		/// <summary>
		/// Helper method for removing a Delegate from one of the vectors in the maps
		/// </summary>
		/// <param name="listeners">Reference to the vector to search and remove from</param>
		/// <param name="delegate">Instance to remove, used to compare to instance in vector</param>
		FORCE_INLINE void RemoveListener(std::vector<Delegate>& listeners, Delegate&& delegate)
		{
			listeners.erase(
				std::remove(listeners.begin(), listeners.end(), std::move(delegate)),
				listeners.end()
			);
		}
	};
}
#endif // !ULTREALITY_UTILITIES_EVENT_DISPATCHER_H
