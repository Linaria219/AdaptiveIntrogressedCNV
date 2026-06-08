#!/bin/bash


if [ $# -ne 4 ]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: $0 <filename> <line_number> <start_column> <end_column>"
    echo "Example: $0 data.txt 5 3 8"
    exit 1
fi

FILE=$1
LINE_NUM=$2
START_COL=$3
END_COL=$4

if [ ! -f "$FILE" ]; then
    echo "Error: File '$FILE' does not exist"
    exit 1
fi

TOTAL_LINES=$(wc -l < "$FILE")
if [ "$LINE_NUM" -lt 1 ] || [ "$LINE_NUM" -gt "$TOTAL_LINES" ]; then
    echo "Error: Line number $LINE_NUM is invalid, file has $TOTAL_LINES lines"
    exit 1
fi

if [ "$START_COL" -lt 1 ] || [ "$END_COL" -lt 1 ] || [ "$START_COL" -gt "$END_COL" ]; then
    echo "Error: Column range $START_COL-$END_COL is invalid"
    exit 1
fi

awk -v line_num="$LINE_NUM" -v start_col="$START_COL" -v end_col="$END_COL" '
NR == line_num {
    if (end_col > NF) {
        print "Warning: End column", end_col, "is beyond the actual column count", NF
        end_col = NF
    }
    if (start_col > NF) {
        print "Error: Start column", start_col, "is beyond the actual column count", NF
        exit 1
    }
    
    min = $(start_col)
    max = $(start_col)
    
    for (i = start_col; i <= end_col; i++) {
        if ($i + 0 == $i) {  # Check if the value is a number
            if ($i < min) min = $i
            if ($i > max) max = $i
        } else {
            print "Warning: Column", i, "is not a number:", $i
        }
    }
    
    printf "File: %s\n", FILENAME
    printf "Line Number: %d\n", line_num
    printf "Column Range: %d-%d\n", start_col, end_col
    printf "Maximum Value: %s\n", max
    printf "Minimum Value: %s\n", min
    printf "Number of Values Processed: %d\n", (end_col - start_col + 1)
    exit  
}

END {
    if (NR < line_num) {
        print "Error: File has only", NR, "lines, cannot read line", line_num
        exit 1
    }
}
' "$FILE"
