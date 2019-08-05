#include <stdio.h>
#include <stdlib.h>
#include "JNIinterface.h"
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <poll.h>

#include "fpga_pci.h"
#include "fpga_mgmt.h"
#include "fpga_dma.h"
#include "utils/lcd.h"

#include "test_dram_dma_common.h"


#define	MEM_16G              (1ULL << 34)
#define USER_INTERRUPTS_MAX  (16)

/* use the standard out logger */
static const struct logger *logger = &logger_stdout;

int buffer_size = 2048 /*Modify*/

char* received;
char* output = malloc(buffer_size);;

JNIEXPORT void JNICALL Java_JNIinterface_JNI_1send_1string_1to_1c(JNIEnv *env, jobject obj, jstring a){
	const char *str = (*env)->GetStringUTFChars(env, a, 0);
	/* printf("Sent \n\n%s\n\nfrom Java to c\n",str);*/
	dma_send(0,buffer_size,str);




	(*env)->ReleaseStringUTFChars(env, a, str);
	return;
}

JNIEXPORT jstring JNICALL Java_JNIinterface_JNI_1request_1from_1c(JNIEnv *env, jobject obj){
	dma_read(0,buffer_size);
	char* out_str = output;
	return (*env)->NewStringUTF(env,out_str);
}


/**
 * This function send data to DDR DIMMS
 * This function fills a buffer with random data and then uses DMA to copy that
 * buffer into each of the 4 DDR DIMMS.
 */
int dma_send(int slot_id, size_t buffer_size, char* data) {
    int write_fd, read_fd, dimm, rc;

    write_fd = -1;
    read_fd = -1;

    uint8_t *write_buffer = malloc(buffer_size);
    if (write_buffer == NULL) {
        rc = -ENOMEM;
        goto out;
    }

    write_fd = fpga_dma_open_queue(FPGA_DMA_XDMA, slot_id,
        /*channel*/ 0, /*is_read*/ false);
    fail_on((rc = (write_fd < 0) ? -1 : 0), out, "unable to open write dma queue");

    rc = fill_buffer(write_buffer, buffer_size, data); // Need to check with Amazon [NO access]
    fail_on(rc, out, "unabled to initialize buffer");

    for (dimm = 0; dimm < 4; dimm++) {
        start = clock();//added by Haoxuan

        rc = fpga_dma_burst_write(write_fd, write_buffer, buffer_size,
            dimm * MEM_16G);

        fail_on(rc, out, "DMA write failed on DIMM: %d", dimm);
    }
		out:
		    if (write_buffer != NULL) {
		        free(write_buffer);
		    }
		    if (read_buffer != NULL) {
		        free(read_buffer);
		    }
		    if (write_fd >= 0) {
		        close(write_fd);
		    }
		    if (read_fd >= 0) {
		        close(read_fd);
		    }
		    /* if there is an error code, exit with status 1 */
		    return (rc != 0 ? 1 : 0);
}

int dma_read(int slot_id, size_t buffer_size) {
    int read_fd, dimm, rc;
    read_fd = -1;
    uint8_t *read_buffer = malloc(buffer_size);
    if (read_buffer == NULL) {
        rc = -ENOMEM;
        goto out;
    }

    read_fd = fpga_dma_open_queue(FPGA_DMA_XDMA, slot_id,
        /*channel*/ 0, /*is_read*/ true);
    fail_on((rc = (read_fd < 0) ? -1 : 0), out, "unable to open read dma queue");

    for (dimm = 0; dimm < 4; dimm++) {
        rc = fpga_dma_burst_read(read_fd, read_buffer, buffer_size,
            dimm * MEM_16G);
        fail_on(rc, out, "DMA read failed on DIMM: %d", dimm);
    }
		sprintf(output, "%d", *read_buffer);

		out:
		    if (write_buffer != NULL) {
		        free(write_buffer);
		    }
		    if (read_buffer != NULL) {
		        free(read_buffer);
		    }
		    if (write_fd >= 0) {
		        close(write_fd);
		    }
		    if (read_fd >= 0) {
		        close(read_fd);
		    }
		    /* if there is an error code, exit with status 1 */
		    return (rc != 0 ? 1 : 0);
}
