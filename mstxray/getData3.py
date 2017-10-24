import sys
import numpy as np
import MDSplus

def getData3(cardNumber, channel):
    HEADER_LENGTH = 3
    FILE_HEADER = 16

    # Read in the .dat file
    try:
        data = open(str(cardNumber) + '.dat', 'r');
    except:
        print "{}.dat not found".format(cardNumber)
        sys.exit(1)

    # Set the channel to an int
    channel = int(channel)
    # Store every 2 bytes (16 bits) into an element of the array 'f'
    f = np.fromfile(data, dtype=np.int16)
    # Initialize a list to hold the data points of the desired channel
    data_points = []
    # "start" is the starting point of the current set of data points being evaluated
    start = 10
    # "offset" is the starting point of the header associated with the current set of
    # data points being evaluated.
    offset = 0

    # Iterate through the entire set of data points and headers, up to the desired channel
    for current_channel in range(1, channel + 1):
        print 'CHANNEL {}: '.format(current_channel),
        print hex(f[offset] & 0xffff)[0:-1] + ', ' + hex(f[offset + 1] & 0xffff)[0:-1]
        if len(f) < (offset + 11):
            print "1: There are not {} channels, try again".format(channel)
            sys.exit()
        # Grab the lower 2 bytes of the number indicating how many samples there are for the channel
        lowerBits = f[offset + 8]
        # Grab the upper 2 bytes of this number
        upperBits = f[offset + 9]
        # Create the number
        sample_length = 2 * (((upperBits & 0xfff) << 16) + lowerBits)
        print "Event Length = {}".format(sample_length)
        # check if the desired channel has been reached
        if(channel == current_channel):
	    # Check if valid channel, if not then terminate program
            if len(f) < (sample_length + start):
                print "2: There are not {} channels, try again".format(channel)
                sys.exit()
            # if so, then store the data points of this channel into
            # the list "data_points"
            for i in range(start, sample_length + start):
                data_points.append(f[i])
            break
 
        # Update the starting point of the data and the header
        start += (sample_length + 10)
        offset += (start - 10)

    return data_points
