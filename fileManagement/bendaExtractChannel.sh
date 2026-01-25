#!/bin/sh
# Extract a channel from an 8 or 16 channel WAV file from TeeGrid.
# Relies on ffmpeg
# Usage: bendaExtractChannel.sh filename.wav channelNum 
# channelNum starts at 0 and goes to either 7 or 15.

# Get command-line arguments. We aren't doing much checking because we
# have faith that if you are using this, you know what you are doing.

fullfile=$1
channelNum=$2
chaNum='c'$2
baseFileName=`basename $fullfile '.wav'`


# We do check to see if you've made a typographical error in the source filename.
if [ ! -f "$fullfile" ]; then
  # Input file does not exist
  echo "Source file $fullfile does not exist. "
  exit 

else

  # We use ffmpeg to extract the channel.  https://www.ffmpeg.org/download.html

  ffmpeg -i $fullfile -af "pan=mono|c0=$chaNum" $baseFileName-$chaNum.wav

  if [ ! -f "$baseFileName-$chaNum.wav" ]; then
    echo "Something went wrong with generating the new file."
    else
    echo "Channel $channelNum successfully extracted to $baseFileName-$chaNum.wav."
  fi

fi

