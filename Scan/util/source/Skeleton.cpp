#include <iostream>

#include <getopt.h>
#include <stdlib.h>

#include "XiaData.hpp"

// Local files
#include "Skeleton.hpp"

// Define the name of the program.
#ifndef PROG_NAME
#define PROG_NAME "Spooky"
#endif

///////////////////////////////////////////////////////////////////////////////
// class skeletonUnpacker
///////////////////////////////////////////////////////////////////////////////

/** Process all events in the event list.
  * \param[in]  addr_ Pointer to a location in memory. 
  * \return Nothing.
  */
void skeletonUnpacker::ProcessRawEvent(ScanInterface *addr_/*=NULL*/){
	if(!addr_){ return; }
	
	XiaData *current_event = NULL;
	
	// Fill the processor event deques with events
	while(!rawEvent.empty()){
		current_event = rawEvent.front();
		rawEvent.pop_front(); // Remove this event from the raw event deque.
		
		// Check that this channel event exists.
		if(!current_event){ continue; }

		// Send the event to the scan interface object for processing.
		if(addr_->AddEvent(current_event))
			addr_->ProcessEvents();
	}
	
	// Finish up with this raw event.
	addr_->ProcessEvents();
}

///////////////////////////////////////////////////////////////////////////////
// class skeletonScanner
///////////////////////////////////////////////////////////////////////////////

/// Default constructor.
skeletonScanner::skeletonScanner() : ScanInterface() {
	init = false;
}

/// Destructor.
skeletonScanner::~skeletonScanner(){
	if(init){
		// Handle some cleanup.
	}
}

/** ExtraCommands is used to send command strings to classes derived
  * from ScanInterface. If ScanInterface receives an unrecognized
  * command from the user, it will pass it on to the derived class.
  * \param[in]  cmd_ The command to interpret.
  * \param[out] arg_ Vector or arguments to the user command.
  * \return True if the command was recognized and false otherwise.
  */
bool skeletonScanner::ExtraCommands(const std::string &cmd_, std::vector<std::string> &args_){
	if(cmd_ == "mycmd"){
		if(args_.size() >= 1){
			// Handle the command.
		}
		else{
			std::cout << msgHeader << "Invalid number of parameters to 'mycmd'\n";
			std::cout << msgHeader << " -SYNTAX- mycmd <param>\n";
		}
	}
	else{ return false; } // Unrecognized command.

	return true;
}

/** ExtraArguments is used to send command line arguments to classes derived
 * from ScanInterface. It has its own instance of getopts to look for its known 
 * parameters. 
 * \param[in] argc : The number of command line arguments
 * \param[in] argv[] : The arrary containing all command line arguments */
void skeletonScanner::ExtraArguments(int argc, char *argv[]) {
    struct option opts[] = {
        { "xtra", no_argument, NULL, 'x'},
        { "kick", required_argument, NULL, 'k'},
        { "ugh", no_argument, NULL, 'u'},
        { NULL, no_argument, NULL, 0 }
    };
    
    std::string optstr = "xk:u";
    int idx = 0;
    int retval = 0;
    
    while ( (retval = getopt_long(argc,argv,optstr.c_str(),opts,&idx)) != -1) {
        switch(retval) {
        case 'x':
            std::cout << "Got option X" << std::endl;
            break;
        case 'k':
            std::cout << "Got option k  " << optarg << std::endl;
            break;
        case ':' :
        case '?' :
        default:
            ArgHelp();
            exit(0);
        }
    }
}

/** CmdHelp is used to allow a derived class to print a help statement about
  * its own commands. This method is called whenever the user enters 'help'
  * or 'h' into the interactive terminal (if available).
  * \param[in]  prefix_ String to append at the start of any output.
  * \return Nothing.
  */
void skeletonScanner::CmdHelp(){
	std::cout << "   mycmd <param> - Do something useful.\n";
}

/** ArgHelp is used to allow a derived class to print a help statment about
  * its own command line arguments. This method is called at the end of
  * the ScanInterface::help method.
  * \return Nothing.
  */
void skeletonScanner::ArgHelp(){
	std::cout << "   --myarg - A useful command line argument.\n";
}

/** SyntaxStr is used to print a linux style usage message to the screen.
  * \param[in]  name_ The name of the program.
  * \return Nothing.
  */
void skeletonScanner::SyntaxStr(char *name_){ 
	std::cout << " usage: " << std::string(name_) << " [input] [options]\n"; 
}

/** Initialize the map file, the config file, the processor handler, 
  * and add all of the required processors.
  * \param[in]  prefix_ String to append to the beginning of system output.
  * \return True upon successfully initializing and false otherwise.
  */
bool skeletonScanner::Initialize(std::string prefix_){
	if(init){ return false; }

	// Do some initialization.

	return (init = true);
}

/** Peform any last minute initialization before processing data. 
  * /return Nothing.
  */
void skeletonScanner::FinalInitialization(){
	// Do some last minute initialization before the run starts.
}

/** Receive various status notifications from the scan.
  * \param[in] code_ The notification code passed from ScanInterface methods.
  * \return Nothing.
  */
void skeletonScanner::Notify(const std::string &code_/*=""*/){
	if(code_ == "START_SCAN"){  }
	else if(code_ == "STOP_SCAN"){  }
	else if(code_ == "SCAN_COMPLETE"){ std::cout << msgHeader << "Scan complete.\n"; }
	else if(code_ == "LOAD_FILE"){ std::cout << msgHeader << "File loaded.\n"; }
	else if(code_ == "REWIND_FILE"){  }
	else{ std::cout << msgHeader << "Unknown notification code '" << code_ << "'!\n"; }
}

/** Return a pointer to the Unpacker object to use for data unpacking.
  * If no object has been initialized, create a new one.
  * \return Pointer to an Unpacker object.
  */
Unpacker *skeletonScanner::GetCore(){ 
	if(!core){ core = (Unpacker*)(new skeletonUnpacker()); }
	return core;
}

/** Add a channel event to the deque of events to send to the processors.
  * This method should only be called from skeletonUnpacker::ProcessRawEvent().
  * \param[in]  event_ The raw XiaData to add to the channel event deque.
  * \return False.
  */
bool skeletonScanner::AddEvent(XiaData *event_){
	if(!event_){ return false; }

	// Handle the individual XiaData. Maybe add it to a detector's event list or something.
	// Do nothing with it for now.
	delete event_;
	
	return false;
}

/** Process all channel events read in from the rawEvent.
  * This method should only be called from skeletonUnpacker::ProcessRawEvent().
  * \return False.
  */
bool skeletonScanner::ProcessEvents(){
	// Process all of the events added so far.
	return false;
}

int main(int argc, char *argv[]){
	// Define a new unpacker object.
	skeletonScanner scanner;
	
	// Set the output message prefix.
	scanner.SetProgramName(std::string(PROG_NAME));	
	
	// Initialize the scanner.
	scanner.Setup(argc, argv);

	// Run the main loop.
	int retval = scanner.Execute();
	
	scanner.Close();
	
	return retval;
}