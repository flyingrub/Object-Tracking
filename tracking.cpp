#include <opencv2/opencv.hpp>
#include <opencv2/tracking.hpp>
#include <opencv2/core/ocl.hpp>
#include <numeric>
using namespace cv;
using namespace std;

// Convert to string
#define SSTR( x ) static_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()


Ptr<Tracker> getTracker(string trackerType) {
    Ptr<Tracker> tracker;

    #if (CV_MINOR_VERSION < 3)
    {
        tracker = Tracker::create(trackerType);
    }
    #else
    {
        if (trackerType == "BOOSTING")
            tracker = TrackerBoosting::create();
        if (trackerType == "MIL")
            tracker = TrackerMIL::create();
        if (trackerType == "KCF")
            tracker = TrackerKCF::create();
        if (trackerType == "TLD")
            tracker = TrackerTLD::create();
        if (trackerType == "MEDIANFLOW")
            tracker = TrackerMedianFlow::create();
        if (trackerType == "GOTURN")
            tracker = TrackerGOTURN::create();
        if (trackerType == "MOSSE")
            tracker = TrackerMOSSE::create();
        if (trackerType == "CSRT")
            tracker = TrackerCSRT::create();
    }
    #endif
    return tracker;
}

int main(int argc, char **argv)
{
    // List of tracker types in OpenCV 3.4.1
    vector<string> trackerTypes = {"BOOSTING", "MIL", "KCF", "TLD","MEDIANFLOW", "MOSSE", "CSRT"};
    // vector <string> trackerTypes(types, std::end(types));
    int currentTrackerType = 0;
    // Create a tracker

    Ptr<Tracker> tracker;
    tracker = getTracker(trackerTypes[currentTrackerType] );
    // Read video
    VideoCapture video("/home/fly/workspace/VR/SequenceAvecFleurs/Image%03d.tif");
    // VideoCapture video("/home/fly/workspace/VR/Ghost4/GITS%03d.bmp");
    Size S = Size((int) video.get(CAP_PROP_FRAME_WIDTH),    // Acquire input size
                  (int) video.get(CAP_PROP_FRAME_HEIGHT));
    int ex = static_cast<int>(video.get(CAP_PROP_FOURCC));     // Get Codec Type- Int form
    VideoWriter outputVideo;                                        // Open the output
    outputVideo.open("test.avi", CV_FOURCC('M','J','P','G'), video.get(CAP_PROP_FPS), S, true);

    // Exit if video is not opened
    if(!video.isOpened())
    {
        cout << "Could not read video file" << endl;
        return 1;
    }

    // Read first frame
    Mat frame;
    bool ok = video.read(frame);

    // Define initial bounding box
    Rect2d originalBox(287, 23, 86, 320);

    // Uncomment the line below to select a different bounding box
    originalBox = selectROI(frame, false);
    Rect2d bbox = originalBox;
    // Display bounding box.
    rectangle(frame, originalBox, Scalar( 255, 0, 0 ), 2, 1 );

    imshow("Tracking", frame);
    tracker->init(frame, originalBox);
    bool hasNext = true;
    vector<float> allFPS;
    while(video.isOpened())
    {
        hasNext = video.read(frame);

        if (!hasNext) {
          double sum = std::accumulate(allFPS.begin(), allFPS.end(), 0.0);
          double mean = sum / allFPS.size();
          cout << trackerTypes[currentTrackerType] << " - meanFPS:" << mean << endl;
          allFPS.clear();

          currentTrackerType++;
          if (currentTrackerType == trackerTypes.size()) {
            break;
          }
          tracker = getTracker(trackerTypes[currentTrackerType]);
          video.set(CV_CAP_PROP_POS_FRAMES, 0);
          hasNext = video.read(frame);
          tracker->init(frame, originalBox);
        }

        // Start timer
        double timer = (double)getTickCount();

        // Update the tracking result
        bool ok = tracker->update(frame, bbox);

        // Calculate Frames per second (FPS)
        float fps = getTickFrequency() / ((double)getTickCount() - timer);
        allFPS.push_back(fps);
        if (ok)
        {
            // Tracking success : Draw the tracked object
            rectangle(frame, bbox, Scalar( 255, 0, 0 ), 2, 1 );
        }
        else
        {
            // Tracking failure detected.
            putText(frame, "Tracking failure detected", Point(100,80), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(0,0,255),2);
        }

        // Display tracker type on frame
        putText(frame, trackerTypes[currentTrackerType]  + " Tracker", Point(100,20), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(50,170,50),2);

        // Display FPS on frame
        putText(frame, "FPS : " + SSTR(int(fps)), Point(100,50), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(50,170,50), 2);

        // Display frame.
        imshow("Tracking", frame);
        outputVideo << frame;

        // Exit if ESC pressed.
        int k = waitKey(1);
        if(k == 27)
        {
            break;
        }
    }
}