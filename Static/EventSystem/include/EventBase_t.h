#ifndef ULTREALITY_UTILITIES_EVENT_BASE_T_H
#define ULTREALITY_UTILITIES_EVENT_BASE_T_H

#include <type_traits>

namespace UltReality::Utilities
{
	/*
		File: EventBase_t.h

		Description:
		Defines a generic base template, EventBase_t, for an event system. This template
		is designed to be the foundation for different event types in an application.
		By templating with an EnumClass type, EventBase_t allows an event system to use
		strongly-types enumerations for various event categories, making event handling
		type-safe and extensible.

		How It Works:
		-	'EventBase_t' is a templated base event structure. It takes an enum class as
			its template parameter ('EnumType'), ensuring type safety by using the EnumClass concept.
		-	Specific event categories can be represented by different enum classes, such as
			'PlatformEvents', which can store OS-related messages like window resize events or
			key presses.
		-	Intermediate classes, such as 'PlatformMessageEvent', can specialize 'EventBase_t'
			with an enum class, providing a common base for related event types.
		-	Concrete event classes, like 'WindowResizeEvent', can inherit from these intermediate
			classes and add fields specific to those events.

		Example Usage:
		-	'PlatformEvents' is defined as an enum class that includes various OS message types
			such as 'WindowResize' or 'KeyPress'. 'PlatformMessageEvent' serves as an intermediate
			class that inherits from 'EventBase_t<PlatformEvents>'.
		-	'WindowResizeEvent' is derived from 'PlatformMessageEvent' and includes additional
			fields for 'width' and 'height', representing a window resize event.
		-	The event handling system can use these classes to process events based on the
			'type' field, dispatching to the appropriate handler based on each event's specific
			enum value.

			// Enum class defining platform message event types
			enum class PlatformEvents
			{
				WindowResize,
				keyPress,
				MouseMove,
				WindowClose,
				...
			};

			// Intermediate class for platform-specific events
			struct PlatformMessageEvent : public EventBase_t<PlatformEvents>
			{
				using EventBase_t::EventBase_t; // Inherit constructors
			};

			// Specialized event class representing a window resize event
			struct WindowResizeEvent : public PlatformMessageEvent
			{
				int width; // New width of the window
				int height; // New height of the window

				// Initializes the event with the specified dimensions, and set the enum type for the event
				WindowResizeEvent(int w, int h) : PlatformMessageEvent(PlatformEvents::WindowResize), width(w), height(h) {}
			};
	*/


	/// <summary>
	/// Concept that verifies a template argument is an enum class type
	/// </summary>
	template<typename T>
	concept CEnumClass = std::is_enum_v<T>;

	/// <summary>
	/// Base template used to create Event types for an event system
	/// </summary>
	/// <typeparam name="EnumType">An enum class type which contains fields representing the type of events supported by the event system</typeparam>
	template<CEnumClass Enum>
	struct EventBase_t
	{
		using EnumType = Enum; // Alias for the enum class type used

		EnumType type; // Stores the type of event this particular instance is, represented by an enum class field from the enum class type provided as a template argument

		explicit EventBase_t(EnumType t) : type(t)
		{}

		/*EventBase_t(const EventBase_t&) = delete;
		EventBase_t& operator=(const EventBase_t&) = delete;
		EventBase_t(EventBase_t&&) = delete;
		EventBase_t& operator=(EventBase_t&&) = delete;*/
	};
}

#endif // !ULTREALITY_UTILITIES_EVENT_BASE_T_H
