/** \file TtuTimingProcessor.cpp
 * \brief Analyzes pulser signals
 *
 * Analyzes pulser signals for electronics and high resolution
 * timing applications.
 *
 * \author S. V. Paulauskas
 * \date 10 July 2009
 */
#include <fstream>
#include <iostream>

#include <cmath>

#include "DammPlotIds.hpp"
#include "TtuTimingProcessor.hpp"
#include "RawEvent.hpp"
#include "TimingCalibrator.hpp"

namespace dammIds {
    namespace ttutiming {
        const int D_TIMEDIFF     = 0;//!< Time Difference
        const int D_PROBLEMSTUFF = 1;//!< Histogram for problems

        const int DD_QDC         = 2;//!< QDCs
        const int DD_MAX         = 3;//!< Max Values
        const int DD_PVSP        = 4;//!< Phase vs. Phase
        const int DD_MAXVSTDIFF  = 5;//!< Maximum Start vs. TDiff
        const int DD_QDCVSMAX    = 6;//!< QDC Start vs. Max Start
        const int DD_AMPMAPSTART = 7;//!< Amplitude Map Start
        const int DD_AMPMAPSTOP  = 8;//!< Amplitude Map Stop
        const int DD_SNRANDSDEV  = 9;//!< SNR and Standard Deviation
        const int DD_PROBLEMS    = 13;//!< 2D problems
        const int DD_MAXSVSTDIFF = 14;//!< Maximum Stop vs. Tdiff
    }
}

using namespace std;
using namespace dammIds::ttutiming;

TtuTimingProcessor::TtuTimingProcessor(): EventProcessor(dammIds::ttutiming::OFFSET, 
    							 dammIds::ttutiming::RANGE, 								 "TtuTimingProcessor") {
    associatedTypes.insert("pulser");
    associatedTypes.insert("vandle");
    associatedTypes.insert("beta");
}

void TtuTimingProcessor::DeclarePlots(void) {
    DeclareHistogram1D(D_TIMEDIFF, SC, "Time Difference");
    DeclareHistogram1D(D_PROBLEMSTUFF, S5, "Problem Stuff");

    DeclareHistogram2D(DD_QDC, SD, S1,"QDC");
    DeclareHistogram2D(DD_MAX, SC, S1, "Max");
    DeclareHistogram2D(DD_PVSP, SC, SC,"Phase vs. Phase");
    //DeclareHistogram2D(DD_MAXVSTDIFF, SC, SC, "Max vs. Time Diff");
    //DeclareHistogram2D(DD_QDCVSMAX, SC, SD,"QDC vs Max");
    //DeclareHistogram2D(DD_AMPMAPSTART, S7, SC,"Amp Map Start");
    //DeclareHistogram2D(DD_AMPMAPSTOP, S7, SC,"Amp Map Stop");
    //DeclareHistogram2D(DD_SNRANDSDEV, S8, S2, "SNR and SDEV R01/L23");
    //DeclareHistogram2D(DD_PROBLEMS, SB, S5, "Problems - 2D");
}

bool TtuTimingProcessor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return false;

    plot(D_PROBLEMSTUFF, 30);

    pulserMap_.clear();

    static const vector<ChanEvent*> & pulserEvents =
        event.GetSummary("pulser")->GetList();

    for(vector<ChanEvent*>::const_iterator itTtuTiming = pulserEvents.begin();
	itTtuTiming != pulserEvents.end(); itTtuTiming++) {
        unsigned int location = (*itTtuTiming)->GetChanID().GetLocation();
        string subType = (*itTtuTiming)->GetChanID().GetSubtype();

        TimingDefs::TimingIdentifier key(location, subType);
        pulserMap_.insert(make_pair(key, HighResTimingData(*(*itTtuTiming))));
    }

    if(pulserMap_.empty() || pulserMap_.size()%2 != 0) {
        plot(D_PROBLEMSTUFF, 27);
	EndProcess();
        return(false);
    }
     
    HighResTimingData start =
        (*pulserMap_.find(make_pair(0,"start"))).second;
    HighResTimingData stop  =
        (*pulserMap_.find(make_pair(0,"stop"))).second;

    // static int counter = 0;
    // for(Trace::const_iterator it = start.GetTrace()->begin();
    //     it!= start.GetTrace()->end(); it++)
    //     plot(DD_PROBLEMS, int(it-start.GetTrace()->begin()), counter, *it);
    // counter ++;


    // cout << "TtuTimingProcessor::Process : "  << start.GetTraceQdc()
    //      << " " << start.GetMaximumPosition()
    //      << " " << start.GetStdDevBaseline()
    //      << " " << start.GetPhase() << endl;
    plot(DD_QDC, start.GetTraceQdc(), 0);
    plot(DD_MAX, start.GetMaximumValue(), 0);

    
    if(start.GetIsValid() && stop.GetIsValid()) {
        double timeDiff = stop.GetHighResTimeInNs() - start.GetHighResTimeInNs();
        double timeRes  = 100.; //500 ps/bin
        double timeOff  = 1000.;
        double phaseX   = 0.;

        // cout << timeDiff*timeRes + timeOff << " "
        //      << start.GetPhase() * timeRes - phaseX << endl;

        // cout << setprecision(13) << "TtuTiming Processor: " << start.GetHighResTimeInNs() << " "
        //      << stop.GetHighResTimeInNs() << " " << timeDiff << " "
        //      << start.GetTrace().GetPhase() << " " << stop.GetTrace().GetPhase() << " " 
	//       << start.GetFilterTime() << " " << stop.GetFilterTime() << " " 
        //      << (double)start.GetPhase()+(double)start.GetFilterTime() << " "
        //      << stop.GetPhase() + stop.GetFilterTime() << " "
        //      << (double)(start.GetPhase()+start.GetFilterTime()) -
        //     (double)(stop.GetPhase() + stop.GetFilterTime())
        //      << endl;

	  plot(D_TIMEDIFF, timeDiff*timeRes + timeOff);
        plot(DD_PVSP, start.GetTrace().GetPhase()*timeRes-phaseX,
            stop.GetTrace().GetPhase()*timeRes-phaseX);

	/*
        plot(DD_MAXVSTDIFF, timeDiff*timeRes+timeOff, start.GetMaximumValue());
        plot(DD_QDCVSMAX, start.GetMaximumValue(), start.GetTraceQdc());
        plot(DD_QDC, stop.GetTraceQdc(), 1);
        plot(DD_MAX, stop.GetMaximumValue(), 1);
        plot(DD_SNRANDSDEV, start.GetTrace().GetSignalToNoiseRatio()+50, 0);
        plot(DD_SNRANDSDEV, start.GetStdDevBaseline()*timeRes+timeOff, 1);
        plot(DD_SNRANDSDEV, stop.GetTrace().GetSignalToNoiseRatio()+50, 2);
        plot(DD_SNRANDSDEV, stop.GetStdDevBaseline()*timeRes+timeOff, 3);
	 */
    }

    EndProcess();
    return(true);
}
