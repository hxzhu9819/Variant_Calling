# Purpsoe
The goal of this part is to convert all required values to binary format. The bitmap is as follows,
```html
Bit map:
Word0:(header 00)             0:1                2:15        16:31        32:47     48:63                    64:288
                              category           --          #reads       #hp       hp_start_addr            addr_of_last_row. [batch size (sram size)]
WordR:(read 01)               0:1                2:15        16:288
                              category           Length      bp+Q (8bit unit)
WordH:(haplotype 10)          0:1                2:15        16:47                            48:79                   80:288
                              category           Length      LOG2_INITIAL_VALUE(32bit)        INITIAL_VALUE(32bit)    bp	
```

## Content
* String.java: uses String to fulfill the goal.
* Untitile.java: NOT IMPLEMENTED

## Log
June 28
* Base line implementation completed via String method (may be slow)
