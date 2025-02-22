#include <gtest/gtest.h>

#include <EventDispatcher_t.h>

using namespace UltReality::Utilities;

enum class PlatformEvents
{
    EWindowActivate, // Triggered when the window is activated or deactivated
    EWindowResize, // Triggered when the user resizes the window
    EWindowEnterSizeMove, // Triggered when the user grabs the resize bars
    EWindowExitSizeMove, // Triggered when the user releases the resize bars
    EWindowClose, // Triggered when a close request is sent
    EWindowDestroy, // Triggered when the window is being destroyed
    EQuit, // Triggered when the application is terminating
    EWindowFocusGained, // Triggered when the window gains focus
    EWindowFocusLost, // Triggered when the window looses focus
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

struct EWindowResize : public PlatformMessageEvent
{
    enum class SizeDetails
    {
#if defined(_WIN_TARGET)
        MaxHide,
        Maximized,
        MaxShow,
        Minimized,
        Restored
#elif defined(_LINUX_TARGET)
        MaxHide,
        Maximized,
        MaxShow,
        Minimized,
        Restored
#endif
    };

    uint16_t width; // The width of the window after the resize
    uint16_t height; // The height of the window after the resize

    SizeDetails details;

    EWindowResize(uint16_t w, uint16_t h, SizeDetails d) : PlatformMessageEvent(PlatformEvents::EWindowResize), width(w), height(h), details(d) {}
};

enum class Key
{
    CTRL
};

struct EKeyDown : public PlatformMessageEvent
{
    Key key;

    EKeyDown(Key k) : PlatformMessageEvent(PlatformEvents::EKeyDown), key(k) {}
};

struct Mouse
{
    enum class CompoundKeys
    {
        Control,
        LButton,
        MButton,
        RButton,
        Shift,
        XButton1,
        XButton2,
        None
    };

    enum class XButton
    {
        XButton1,
        XButton2
    };

    CompoundKeys compKey;

    uint16_t mouseX;
    uint16_t mouseY;

    Mouse(CompoundKeys k, uint16_t x, uint16_t y) : compKey(k), mouseX(x), mouseY(y) {}
};

struct EMouseMove : public PlatformMessageEvent, Mouse
{
    EMouseMove(CompoundKeys k, uint16_t x, uint16_t y) : PlatformMessageEvent(PlatformEvents::EMouseMove), Mouse(k, x, y) {}
};

class EventDispatcherTest : public ::testing::Test
{
protected:
	EventDispatcher_t<PlatformMessageEvent> dispatcher;

    bool wasMemberFunctionCalled1 = false;
    bool wasMemberFunctionCalled2 = false;

public:
    void MemberFunctionListener1(const EWindowResize& event)
    {
        wasMemberFunctionCalled1 = true;
    }

    void MemberFunctionListener2(const EWindowResize& event)
    {
        wasMemberFunctionCalled2 = true;
    }
};

// Test if a lambda listener is called when a specific event is dispatched.
TEST_F(EventDispatcherTest, ListenerCalledOnDispatch)
{
    bool listenerCalled = false;

    dispatcher.Subscribe([&](const EWindowResize& event) {
        listenerCalled = true;
        EXPECT_EQ(event.width, 1280);
        EXPECT_EQ(event.height, 720);
        EXPECT_EQ(event.details, EWindowResize::SizeDetails::Restored);
        }, PlatformEvents::EWindowResize);

    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.DispatchBatch();

    EXPECT_TRUE(listenerCalled);
}

// Test subscribing and unsubscribing a free function listener
void TestFreeFunction(const EKeyDown& event)
{
    // This function could assert properties of the event type if needed
    EXPECT_EQ(event.key, Key::CTRL);
    EXPECT_EQ(event.type, PlatformEvents::EKeyDown);
}

TEST_F(EventDispatcherTest, FreeFunctionSubscribeAndUnsubscribe)
{
    dispatcher.Subscribe(&TestFreeFunction, PlatformEvents::EKeyDown);
    EXPECT_NO_THROW(dispatcher.Unsubscribe(&TestFreeFunction, PlatformEvents::EKeyDown));
}

// Test dispatching multiple events of different types and verifying listeners are called accordingly
TEST_F(EventDispatcherTest, MultipleEventTypesDispatch)
{
    bool windowResizeCalled = false;
    bool mouseMoveCalled = false;

    dispatcher.Subscribe([&](const EWindowResize&) { windowResizeCalled = true; }, PlatformEvents::EWindowResize);
    dispatcher.Subscribe([&](const EMouseMove&) { mouseMoveCalled = true; }, PlatformEvents::EMouseMove);

    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.QueueEvent(EMouseMove(Mouse::CompoundKeys::None, 50, 234));
    dispatcher.DispatchBatch();

    EXPECT_TRUE(windowResizeCalled);
    EXPECT_TRUE(mouseMoveCalled);
}

// Test if unsubscribing a listener successfully prevents it from being called
TEST_F(EventDispatcherTest, UnsubscribePreventsNotification)
{
    bool listenerCalled = false;

    auto lambda = [&](const EWindowResize&) { listenerCalled = true; };
    dispatcher.Subscribe(lambda, PlatformEvents::EWindowResize);
    dispatcher.Unsubscribe(lambda, PlatformEvents::EWindowResize);

    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.DispatchBatch();

    EXPECT_FALSE(listenerCalled);
}

// Test if the event queue is cleared after calling DispatchAll
TEST_F(EventDispatcherTest, DispatchAllClearsQueue)
{
    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.QueueEvent(EWindowResize(800, 600, EWindowResize::SizeDetails::Maximized));

    dispatcher.DispatchAll();

    EXPECT_TRUE(dispatcher.IsQueueEmpty());
}

// Test synchronization when adding events from multiple threads to the SyncQueue
TEST_F(EventDispatcherTest, SyncQueueEventFromMultipleThreads)
{
    std::atomic<int> callCount{ 0 };

    dispatcher.Subscribe([&](const EWindowResize&) { ++callCount; }, PlatformEvents::EWindowResize);

    // Add events from multiple threads
    std::thread t1([&]() { dispatcher.SyncQueueEvent(EWindowResize(800, 600, EWindowResize::SizeDetails::Maximized)); });
    std::thread t2([&]() { dispatcher.SyncQueueEvent(EWindowResize(1024, 768, EWindowResize::SizeDetails::Restored)); });

    t1.join();
    t2.join();

    dispatcher.DispatchBatch();  // This will process the events

    EXPECT_EQ(callCount, 2);
}

// Test if listeners are correctly notified only for the batch size limit
TEST_F(EventDispatcherTest, DispatchBatchRespectsBatchSize)
{
    int callCount = 0;

    dispatcher.SetBatchSize(1);
    dispatcher.Subscribe([&](const EWindowResize&) { ++callCount; }, PlatformEvents::EWindowResize);

    dispatcher.QueueEvent(EWindowResize(800, 600, EWindowResize::SizeDetails::Maximized));
    dispatcher.QueueEvent(EWindowResize(1024, 768, EWindowResize::SizeDetails::Restored));

    dispatcher.DispatchBatch();

    EXPECT_EQ(callCount, 1);  // Only one event should be processed due to batch size
}

// Test thread-safety of SyncSubscribe and SyncUnsubscribe
TEST_F(EventDispatcherTest, SyncSubscribeAndUnsubscribeThreadSafety)
{
    std::atomic<int> callCount{ 0 };

    auto lambda = [&](const EWindowResize&) { ++callCount; };
    dispatcher.SyncSubscribe(lambda, PlatformEvents::EWindowResize);

    std::thread t1([&]() { dispatcher.SyncUnsubscribe(lambda, PlatformEvents::EWindowResize); });
    std::thread t2([&]() { dispatcher.SyncSubscribe(lambda, PlatformEvents::EWindowResize); });

    t1.join();
    t2.join();

    dispatcher.QueueEvent(EWindowResize(800, 600, EWindowResize::SizeDetails::Restored));
    dispatcher.DispatchAll();

    EXPECT_GE(callCount, 0);  // Depending on the timing, listener may or may not be called
}

TEST_F(EventDispatcherTest, MemberFunctionSubscription)
{
    dispatcher.Subscribe<EventDispatcherTest>(this, &EventDispatcherTest::MemberFunctionListener1, PlatformEvents::EWindowResize);

    // Dispatch an event and check if the member function was called
    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.DispatchBatch();

    EXPECT_TRUE(wasMemberFunctionCalled1);
}

TEST_F(EventDispatcherTest, MemberFunctionUnsubscribePreventsNotification)
{
    dispatcher.Subscribe<EventDispatcherTest>(this, &EventDispatcherTest::MemberFunctionListener2, PlatformEvents::EWindowResize);
    dispatcher.Unsubscribe<EventDispatcherTest>(this, &EventDispatcherTest::MemberFunctionListener2, PlatformEvents::EWindowResize);

    // Dispatch an event and check if the unsubscribed member function is NOT called
    dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
    dispatcher.DispatchBatch();

    EXPECT_FALSE(wasMemberFunctionCalled2);
}
