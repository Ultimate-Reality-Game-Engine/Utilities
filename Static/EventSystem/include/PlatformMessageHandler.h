#ifndef ULTREALITY_UTILITIES_PLATFORM_MSG_HANDLER_H
#define ULTREALITY_UTILITIES_PLATFORM_MSG_HANDLER_H

#include <EventBase_t.h>
#include <EventDispatcher.h>

#if defined(_WIN_TARGET)
#include <windows.h>
#elif defined(_LINUX_TARGET)

#elif defined(_MAC_TARGET)

#endif

namespace UltReality::Utilities
{
	/// <summary>
	/// Enum class defining the events supported at the platform level
	/// </summary>
	enum class PlatformEvents
	{
		EWindowActivate, // Triggered when the window is activated or deactivated
		EWindowResize, // Triggered when the user resizes the window
		EWindowEnterSizeMove, // Triggered when the user grabs the resize bars
		EWindowExitSizeMove, // Triggered when the user releases the resize bars
		EWindowDestroy, // Triggered when the window is being destroyed
		EWindowMenuChar, // Triggered when a menu is active and the user presses a key that does not correspond to any mnemonic or accelerator key
		EWindowGetMinMaxInfo, // Triggered to provide size of window
		EMouseLDown, // Triggered when left mouse button pressed
		EMouseLUp, // Triggered when left mouse button released
		EMouseLDClick, // Triggered when left mouse button is double clicked
		EMouseMDown, // Triggered when middle mouse button pressed
		EMouseMUp, // Triggered when middle mouse button released
		EMouseMDClick, // Triggered when middle mouse button double clicked
		EMouseRDown, // Triggered when right mouse button pressed
		EMouseRUp, // Triggered when right mouse button released
		EMouseRDClick, // Triggered when right mouse button double clicked
		EMouseXDClick, // Triggered when one of the X buttons on the mouse is double clicked
		EMouseXDown, // Triggered when one of the X buttons on the mouse is pressed
		EMouseXUp, // Triggered when one of the X buttons on the mouse is released
		EMouseHover, // Triggered when the cursor hovers over the client area of the window
		EMouseLeave, // Triggered when the cursor leaves the client area of the window
		EMouseMove, // Triggered when the mouse moved
		EMouseWheel, // Triggered when the mouse wheel is rotated
		EKeyDown, // Triggered when keyboard key pressed
		EKeyUp, // Triggered when keyboard key released
		ESysKeyDown, // Triggered when user presses F10 or holds down the ALT key and presses another key
		ESysKeyUp // Triggered when the user releases a key that was pressed while the ALT key was held down
	};

	/// <summary>
	/// Base event for the platform level
	/// </summary>
	struct PlatformMessageEvent : public EventBase_t<PlatformEvents>
	{
		using EventBase_t::EventBase_t; // Inherit constructors
	};

	struct EWindowActivate : public PlatformMessageEvent
	{
		enum class ActivationDetails
		{
#if defined(_WIN_TARGET)
			Active = WA_ACTIVE, // Activated by some method other than a mouse click
			ClickActive = WA_CLICKACTIVE, // Activated by a mouse click
			Inactive = WA_INACTIVE // Deactivated
#endif
		};

		ActivationDetails details;

		EWindowActivate(ActivationDetails d) : PlatformMessageEvent(PlatformEvents::EWindowActivate), details(d) {}
	};

	struct EWindowResize : public PlatformMessageEvent
	{
		enum class SizeDetails
		{
#if defined(_WIN_TARGET)
			MaxHide = SIZE_MAXHIDE, // Set when resize event indicates that some other window was just maximized
			Maximized = SIZE_MAXIMIZED, // Set when resize event was to maximize the window
			MaxShow = SIZE_MAXSHOW, // Set when resize event indicates that some other window was just restored to its former size
			Minimized = SIZE_MINIMIZED, // Set when resize event was to minimize the window
			Restored = SIZE_RESTORED // Set when resize event was to restore minimized window
#endif
		};

		uint16_t width; // The width of the window after the resize
		uint16_t height; // The height of the window after the resize

		SizeDetails details;

		EWindowResize(uint16_t w, uint16_t h, SizeDetails d) : PlatformMessageEvent(PlatformEvents::EWindowResize), width(w), height(h), details(d) {}
	};

	struct EWindowEnterSizeMove : public PlatformMessageEvent
	{
		EWindowEnterSizeMove() : PlatformMessageEvent(PlatformEvents::EWindowEnterSizeMove) {}
	};

	struct EWindowExitSizeMove : public PlatformMessageEvent
	{
		EWindowExitSizeMove() : PlatformMessageEvent(PlatformEvents::EWindowExitSizeMove) {}
	};

	struct EWindowDestroy : public PlatformMessageEvent
	{
		EWindowDestroy() : PlatformMessageEvent(PlatformEvents::EWindowDestroy) {}
	};

	struct EWindowMenuChar : public PlatformMessageEvent
	{
		EWindowMenuChar() : PlatformMessageEvent(PlatformEvents::EWindowMenuChar) {}
	};

	struct EWindowMinMax : public PlatformMessageEvent
	{
		using PMINMAXINFO =
#if defined(_WIN_TARGET)
			MINMAXINFO*
#elif defined(_LINUX_TARGET)

#elif defined(_MAC_TARGET)

#endif
			;

		PMINMAXINFO info;

		EWindowMinMax(PMINMAXINFO i) : PlatformMessageEvent(PlatformEvents::EWindowGetMinMaxInfo), info(i) {}
	};

	struct Mouse
	{
		enum class CompoundKeys
		{
#if defined(_WIN_TARGET)
			Control = MK_CONTROL, // The CTRL key is down
			LButton = MK_LBUTTON, // The left mouse button is down
			MButton = MK_MBUTTON, // The middle mouse button is down
			RButton = MK_RBUTTON, // The right mouse button is down
			Shift = MK_SHIFT, // The shift key is down
			XButton1 = MK_XBUTTON1, // The first X button is down
			XButton2 = MK_XBUTTON2, // The second X button is down
#endif
			None
		};

		enum class XButton
		{
#if defined(_WIN_TARGET)
			XButton1 = XBUTTON1, // The first X button on the mouse
			XButton2 = XBUTTON2 // The second X button on the mouse
#endif
		};

		CompoundKeys compKey;

		uint16_t mouseX;
		uint16_t mouseY;

		Mouse(CompoundKeys k, uint16_t x, uint16_t y) : compKey(k), mouseX(x), mouseY(y) {}
	};

	struct EMouseLDown : public PlatformMessageEvent, Mouse
	{
		EMouseLDown(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseLDown), Mouse(k, x, y) {}
	};

	struct EMouseLUp : public PlatformMessageEvent, Mouse
	{
		EMouseLUp(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseLUp), Mouse(k, x, y) {}
	};

	struct EMouseLDClick : public PlatformMessageEvent, Mouse
	{
		EMouseLDClick(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseLDClick), Mouse(k, x, y) {}
	};

	struct EMouseMDown : public PlatformMessageEvent, Mouse
	{
		EMouseMDown(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseMDown), Mouse(k, x, y) {}
	};

	struct EMouseMUp : public PlatformMessageEvent, Mouse
	{
		EMouseMUp(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseMUp), Mouse(k, x, y) {}
	};

	struct EMouseMDClick : public PlatformMessageEvent, Mouse
	{
		EMouseMDClick(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseMDClick), Mouse(k, x, y) {}
	};

	struct EMouseRDown : public PlatformMessageEvent, Mouse
	{
		EMouseRDown(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseRDown), Mouse(k, x, y) {}
	};

	struct EMouseRUp : public PlatformMessageEvent, Mouse
	{
		EMouseRUp(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseRUp), Mouse(k, x, y) {}
	};

	struct EMouseRDClick : public PlatformMessageEvent, Mouse
	{
		EMouseRDClick(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseRDClick), Mouse(k, x, y) {}
	};

	struct EMouseXDClick : public PlatformMessageEvent, Mouse
	{
		XButton button;

		EMouseXDClick(XButton b, CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseXDClick), Mouse(k, x, y), button(b) {}
	};

	struct EMouseXDown : public PlatformMessageEvent, Mouse
	{
		XButton button;

		EMouseXDown(XButton b, CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseXDown), Mouse(k, x, y), button(b) {}
	};

	struct EMouseXUp : public PlatformMessageEvent, Mouse
	{
		XButton button;

		EMouseXUp(XButton b, CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseXUp), Mouse(k, x, y), button(b) {}
	};

	struct EMouseHover : public PlatformMessageEvent, Mouse
	{
		EMouseHover(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseHover), Mouse(k, x, y) {}
	};

	struct EMouseLeave : public PlatformMessageEvent
	{
		EMouseLeave() : PlatformMessageEvent(PlatformEvents::EMouseLeave) {}
	};

	struct EMouseMove : public PlatformMessageEvent, Mouse
	{
		EMouseMove(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseMove), Mouse(k, x, y) {}
	};

	struct EMouseWheel : public PlatformMessageEvent, Mouse
	{
#if defined(_WIN_TARGET)
		static constexpr int8_t WheelDeltaMult = 120; // Multiple to be applied to the degrees field. All degrees are in terms of this multiple
#endif
		int16_t degrees;

		EMouseWheel(int16_t d, CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseWheel), Mouse(k, x, y), degrees(d) {}
	};

	enum class Key
	{
#if defined(_WIN_TARGET)
		MouseLButton = VK_LBUTTON, // Left mouse button
		MouseRButton = VK_RBUTTON, // Right mouse button
		Cancel = VK_CANCEL, // Control-break processing
		MouseMButton = VK_MBUTTON, // Middle mouse button
		MouseXButton1 = VK_XBUTTON1, // X1 mouse button
		MouseXButton2 = VK_XBUTTON2, // X2 mouse button
		Back = VK_BACK, // Backspace key
		Tab = VK_TAB, // Tab key
		Clear = VK_CLEAR, // Clear key
		Return = VK_RETURN, // Enter key
		Enter = VK_RETURN, // Enter key
		Shift = VK_SHIFT, // Shift key
		Control = VK_CONTROL, // CTRL key
		CTRL = VK_CONTROL, // CTRL key
		Menu = VK_MENU, // ALT key
		ALT = VK_MENU, // ALT key
		Pause = VK_PAUSE, // Pause key
		Capital = VK_CAPITAL, // Caps Lock key
		CAPS = VK_CAPITAL, // Caps Lock key
		KANA = VK_KANA, // IME Kana mode
		HANGUL = VK_HANGUL, // IME Hangul mode
		IME_ON = VK_IME_ON, // IME On
		JUNJA = VK_JUNJA, // IME Junja mode
		FINAL = VK_FINAL, // IME Final mode
		HANJA = VK_HANJA, // IME Hanja mode
		KANJI = VK_KANJI, // IME Kanji mode
		IME_OFF = VK_IME_OFF, // IME Off
		Escape = VK_ESCAPE, // ESC key
		ESC = VK_ESCAPE, // ESC key
		CONVERT = VK_CONVERT, // IME Convert
		NONCONVERT = VK_NONCONVERT, // IME Non-convert
		ACCEPT = VK_ACCEPT, // IME Accept
		MODECHANGE = VK_MODECHANGE, // IME mode change request
		Space = VK_SPACE, // Spacebar key
		Prior = VK_PRIOR, // Page Up key
		PAGE_UP = VK_PRIOR, // Page Up key
		Next = VK_NEXT, // Page Down key
		PAGE_DOWN = VK_NEXT, // Page Down key
		End = VK_END, // End key
		Home = VK_HOME, // Home key
		Left = VK_LEFT, // Left arrow key
		Up = VK_UP, // Up arrow key
		Right = VK_RIGHT, // Right arrow key
		Down = VK_DOWN, // Down arrow key
		Select = VK_SELECT, // Select key
		Print = VK_PRINT, // Print key
		Execute = VK_EXECUTE, // Execute key
		Snapshot = VK_SNAPSHOT, // Print Screen key
		PRINT_SCREEN = VK_SNAPSHOT, // Print Screen key
		Insert = VK_INSERT, // INS key
		INS = VK_INSERT, // INS key
		Delete = VK_DELETE, // Delete key
		DEL = VK_DELETE, // Delete key
		Help = VK_HELP, // Help key
		Zero = 0x30, // 0 key
		D0 = 0x31, // 1 key
		D2 = 0x32, // 2 key
		D3 = 0x33, // 3 key
		D4 = 0x34, // 4 key
		D5 = 0x35, // 5 key
		D6 = 0x36, // 6 key
		D7 = 0x37, // 7 key
		D8 = 0x38, // 8 key
		D9 = 0x39, // 9 key
		A = 0x41, // A key
		B = 0x42, // B key
		C = 0x43, // C key
		D = 0x44, // D key
		E = 0x45, // E key
		F = 0x46, // F key
		G = 0x47, // G key
		H = 0x48, // H key
		I = 0x49, // I key
		J = 0x4A, // J key
		K = 0x4B, // K key
		L = 0x4C, // L key
		M = 0x4D, // M key
		N = 0x4E, // N key
		O = 0x4F, // O key
		P = 0x50, // P key
		Q = 0x51, // Q key
		R = 0x52, // R key
		S = 0x53, // S key
		T = 0x54, // T key
		U = 0x55, // U key
		V = 0x56, // V key
		W = 0x57, // W key
		X = 0x58, // X key
		Y = 0x59, // Y key
		Z = 0x5A, // Z key
		LWIN = VK_LWIN, // Left Windows key
		RWIN = VK_RWIN, // Right Windows key
		APPS = VK_APPS, // Applications key
		Sleep = VK_SLEEP, // Computer sleep key
		NumPad0 = VK_NUMPAD0, // Numeric keypad 0
		NumPad1 = VK_NUMPAD1, // Numeric keypad 1
		NumPad2 = VK_NUMPAD2, // Numeric keypad 2
		NumPad3 = VK_NUMPAD3, // Numeric keypad 3
		NumPad4 = VK_NUMPAD4, // Numeric keypad 4
		NumPad5 = VK_NUMPAD5, // Numeric keypad 5
		NumPad6 = VK_NUMPAD6, // Numeric keypad 6
		NumPad7 = VK_NUMPAD7, // Numeric keypad 7
		NumPad8 = VK_NUMPAD8, // Numeric keypad 8
		NumPad9 = VK_NUMPAD9, // Numeric keypad 9
		Multiply = VK_MULTIPLY, // Multiply key
		Add = VK_ADD, // Add key
		Separator = VK_SEPARATOR, // Separator key
		Subtract = VK_SUBTRACT, // Subtract key
		Decimal = VK_DECIMAL, // Decimal key
		Divide = VK_DIVIDE, // Divide key
		F1 = VK_F1, // F1 key
		F2 = VK_F2, // F2 key
		F3 = VK_F3, // F3 key
		F4 = VK_F4, // F4 key
		F5 = VK_F5, // F5 key
		F6 = VK_F6, // F6 key
		F7 = VK_F7, // F7 key
		F8 = VK_F8, // F8 key
		F9 = VK_F9, // F9 key
		F10 = VK_F10, // F10 key
		F11 = VK_F11, // F11 key
		F12 = VK_F12, // F12 key
		F13 = VK_F13, // F13 key
		F14 = VK_F14, // F14 key
		F15 = VK_F15, // F15 key
		F16 = VK_F16, // F16 key
		F17 = VK_F17, // F17 key
		F18 = VK_F18, // F18 key
		F19 = VK_F19, // F19 key
		F20 = VK_F20, // F20 key
		F21 = VK_F21, // F21 key
		F22 = VK_F22, // F22 key
		F23 = VK_F23, // F23 key
		F24 = VK_F24, // F24 key
		Numlock = VK_NUMLOCK, // Numlock key
		SCROLL_LOCK = VK_SCROLL, // Scroll lock key
		LShift = VK_LSHIFT, // Left shift key
		RShift = VK_RSHIFT, // Right shift key
		LControl = VK_LCONTROL, // Left CTRL key
		LCTRL = VK_LCONTROL, // Left CTRL key
		RControl = VK_RCONTROL, // Right CTRL key
		RCTRL = VK_RCONTROL, // Right CTRL key
		LMenu = VK_LMENU, // Left ALT key
		LALT = VK_LMENU, // Left ALT key
		RMenu = VK_RMENU, // Right ALT key
		RALT = VK_RMENU, // Right ALT key
		BROWSER_BACK = VK_BROWSER_BACK, // Browser back button
		BROWSER_FORWARD = VK_BROWSER_FORWARD, // Browser forward key
		BROWSER_REFRESH = VK_BROWSER_REFRESH, // Browser refresh key
		BROWSER_STOP = VK_BROWSER_STOP, // Browser stop key
		BROWSER_SEARCH = VK_BROWSER_SEARCH, // Browser search key
		BROWSER_FAVORITES = VK_BROWSER_FAVORITES, // Browser favorites key
		BROWSER_HOME = VK_BROWSER_HOME, // Browser home key
		VolMute = VK_VOLUME_MUTE, // Volume mute key
		VolDown = VK_VOLUME_DOWN, // Volume down key
		VolUp = VK_VOLUME_UP, // Volume up key
#endif
	};

	struct EKeyDown : public PlatformMessageEvent
	{
		Key key;

		EKeyDown(Key k) : PlatformMessageEvent(PlatformEvents::EKeyDown), key(k) {}
	};

	struct EKeyUp : public PlatformMessageEvent
	{
		Key key;

		EKeyUp(Key k) : PlatformMessageEvent(PlatformEvents::EKeyUp), key(k) {}
	};

	struct ESysKeyDown : public PlatformMessageEvent
	{
		Key key;

		ESysKeyDown(Key k) : PlatformMessageEvent(PlatformEvents::ESysKeyDown), key(k) {}
	};

	struct ESysKeyUp : public PlatformMessageEvent
	{
		Key key;

		ESysKeyUp(Key k) : PlatformMessageEvent(PlatformEvents::ESysKeyUp), key(k) {}
	};

	/// <summary>
	/// Class the provides a MsgProc method to be used as the platform's message handler
	/// </summary>
	class PlatformMessageHandler
	{
	private:
		EventDispatcher<PlatformMessageEvent>& m_eventDispatcher;

	public:
		explicit PlatformMessageHandler(EventDispatcher<PlatformMessageEvent>& eventDispatcher) noexcept;

#if defined(_WIN_TARGET)
		/// <summary>
		/// Message handler to be used when the target platform is Windows. Should not be visible otherwise
		/// </summary>
		/// <param name="hwnd">Handle to the window instance</param>
		/// <param name="msg">OS message sent to application</param>
		/// <param name="wParam">First parameter associated with the message</param>
		/// <param name="lParam">Second parameter associated with the message</param>
		/// <returns>Windows defined status or message back to OS</returns>
		LRESULT MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);
#elif defined(_LINUX_TARGET)

#elif defined(_MAC_TARGET)

#endif
	};
}

#endif // !ULTREALITY_UTILITIES_PLATFORM_MSG_HANDLER_H
