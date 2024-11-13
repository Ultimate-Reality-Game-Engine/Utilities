#include <gtest/gtest.h>

#include <EventDispatcher.h>
#include <PlatformMessageHandler.h>

using namespace UltReality::Utilities;

class EventDispatcherTest : public ::testing::Test
{
protected:
	EventDispatcher<PlatformMessageEvent> dispatcher;
};

TEST_F(EventDispatcherTest, ListenerCalledOnDispatch)
{
	bool listenerCalled = false;

	dispatcher.Subscribe<PlatformEvents::EWindowResize>([&](const EWindowResize& event) {
		listenerCalled = true;
		EXPECT_EQ(event.width, 1280);
		EXPECT_EQ(event.width, 720);
		EXPECT_EQ(event.details, EWindowResize::SizeDetails::Restored);
	});

	dispatcher.QueueEvent(EWindowResize(1280, 720, EWindowResize::SizeDetails::Restored));
	dispatcher.DispatchBatch();

	EXPECT_TRUE(listenerCalled);
}