#include <PlatformMessageHandler.h>

namespace UltReality::Utilities
{
	EventDispatcher<PlatformMessageEvent>& PlatformMessageHandler::GetPlatformEventDispatcher() noexcept
	{
		return m_eventDispatcher;
	}

#if defined(_WIN_TARGET)
#include <windowsx.h>

	void PlatformMessageHandler::ProcessPlatformMessages()
	{
		MSG msg = { 0 };

		// Process all available messages in the process queue
		while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE))
		{
			DispatchMessage(&msg);
		}
	}

	void PlatformMessageHandler::PostQuit(const EWindowDestroy&)
	{
		PostQuitMessage(0);
	}

	LRESULT PlatformMessageHandler::MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
	{
		switch (msg)
		{
		// WM_ACTIVATE is sent when the window is activated or deactivated
		case WM_ACTIVATE:
			m_eventDispatcher.QueueEvent(EWindowActivate(static_cast<EWindowActivate::ActivationDetails>(wParam)));

			return 0;

		// WM_SIZE is sent when the user resizes the window
		case WM_SIZE:
			using Details = EWindowResize::SizeDetails;

			m_eventDispatcher.QueueEvent(EWindowResize(LOWORD(lParam), HIWORD(lParam), static_cast<Details>(wParam)));
			
			return 0;

		// WM_ENTERSIZEMOVE is sent when the user grabs the resize bars
		case WM_ENTERSIZEMOVE:
			m_eventDispatcher.QueueEvent(EWindowEnterSizeMove());

			return 0;

		// WM_EXITSIZEMOVE is sent when the user releases the resize bars
		case WM_EXITSIZEMOVE:
			m_eventDispatcher.QueueEvent(EWindowExitSizeMove());

			return 0;

		case WM_CLOSE:
			m_eventDispatcher.QueueEvent(EWindowClose());

			return 0;

		// WM_DESTROY is sent when the window is being destroyed
		case WM_DESTROY:
			m_eventDispatcher.QueueEvent(EWindowDestroy());

			return 0;

		case WM_SETFOCUS:
			m_eventDispatcher.QueueEvent(EWindowFocusGained());

			return 0;

		case WM_KILLFOCUS:
			m_eventDispatcher.QueueEvent(EWindowFocusLost());

			return 0;

		//// The WM_MENUCHAR message is sent when a menu is active and the user presses 
		//// a key that does not correspond to any mnemonic or accelerator key 
		//case WM_MENUCHAR:
		//	m_eventDispatcher.QueueEvent(EWindowMenuChar());

		//	// Don't beep when we alt-enter.
		//	return MAKELRESULT(0, MNC_CLOSE);

		// Catch this message to set maximum and minimum window size available when dragging borders
		case WM_GETMINMAXINFO:
			m_eventDispatcher.QueueEvent(EWindowMinMax((MINMAXINFO*)lParam));

			return 0;

		case WM_LBUTTONDBLCLK:
			using CompKey = Mouse::CompoundKeys;
			m_eventDispatcher.QueueEvent(EMouseLDClick(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_LBUTTONDOWN:
			m_eventDispatcher.QueueEvent(EMouseLDown(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
			
			return 0;

		case WM_LBUTTONUP:
			m_eventDispatcher.QueueEvent(EMouseLUp(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MBUTTONDBLCLK:
			m_eventDispatcher.QueueEvent(EMouseMDClick(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MBUTTONDOWN:
			m_eventDispatcher.QueueEvent(EMouseMDown(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MBUTTONUP:
			m_eventDispatcher.QueueEvent(EMouseMUp(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_RBUTTONDBLCLK:
			m_eventDispatcher.QueueEvent(EMouseRDClick(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_RBUTTONDOWN:
			m_eventDispatcher.QueueEvent(EMouseRDown(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_RBUTTONUP:
			m_eventDispatcher.QueueEvent(EMouseRUp(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MOUSEHOVER:
			m_eventDispatcher.QueueEvent(EMouseHover(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MOUSELEAVE:
			m_eventDispatcher.QueueEvent(EMouseLeave());

			return 0;

		case WM_MOUSEMOVE:
			m_eventDispatcher.QueueEvent(EMouseMove(static_cast<CompKey>(wParam), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_MOUSEWHEEL:
			m_eventDispatcher.QueueEvent(EMouseWheel(GET_WHEEL_DELTA_WPARAM(wParam), static_cast<CompKey>(GET_KEYSTATE_WPARAM(wParam)), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_XBUTTONDBLCLK:
			using XBtn = Mouse::XButton;

			m_eventDispatcher.QueueEvent(EMouseXDClick(static_cast<XBtn>(GET_XBUTTON_WPARAM(wParam)), static_cast<CompKey>(GET_KEYSTATE_WPARAM(wParam)), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_XBUTTONDOWN:
			m_eventDispatcher.QueueEvent(EMouseXDown(static_cast<XBtn>(GET_XBUTTON_WPARAM(wParam)), static_cast<CompKey>(GET_KEYSTATE_WPARAM(wParam)), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_XBUTTONUP:
			m_eventDispatcher.QueueEvent(EMouseXUp(static_cast<XBtn>(GET_XBUTTON_WPARAM(wParam)), static_cast<CompKey>(GET_KEYSTATE_WPARAM(wParam)), GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));

			return 0;

		case WM_KEYDOWN:
			m_eventDispatcher.QueueEvent(EKeyDown(static_cast<Key>(wParam)));

			return 0;

		case WM_KEYUP:
			m_eventDispatcher.QueueEvent(EKeyUp(static_cast<Key>(wParam)));

			return 0;

		case WM_SYSKEYDOWN:
			m_eventDispatcher.QueueEvent(ESysKeyDown(static_cast<Key>(wParam)));

			return 0;

		case WM_SYSKEYUP:
			m_eventDispatcher.QueueEvent(ESysKeyUp(static_cast<Key>(wParam)));

			return 0;
		}

		return DefWindowProc(hwnd, msg, wParam, lParam);
	}
#endif
}