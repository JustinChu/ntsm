/*
 * ProdConKseqRunner.hpp
 *
 *  Created on: Oct 8, 2021
 *      Author: cjustin
 */

#ifndef PRODCONKSEQRUNNER_HPP_
#define PRODCONKSEQRUNNER_HPP_

#include <concurrentqueue.h>
#include <omp.h>
#include <string>
#include <vector>
#include <stdio.h>

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/
#include "kseq_util.h"

using namespace std;

/*
 * Given set of files, object (assume object has process function)
 * Assumes all process functions are thread safe
 */
template<typename T>
void run(const vector<string> &filenames, T &consumer, unsigned threadTotal =
		omp_get_num_threads()) {

	static const unsigned bulkSize = 256;

	assert(threadTotal > filenames.size() + 1);

	const uint64_t bufferSize = threadTotal * bulkSize;
	const uint64_t recycleBufferSize = bufferSize * 2;
	uint64_t readCount;
	uint64_t processedReadsCount;

	unsigned finished = 0;

	moodycamel::ConcurrentQueue<kseq_t> workQueue(bufferSize);
	moodycamel::ConcurrentQueue<kseq_t> recycleQueue(recycleBufferSize);
	typedef std::vector<kseq_t>::iterator iter_t;
	//fill recycleQueue with empty objects
	{
		std::vector<kseq_t> buffer(recycleBufferSize, kseq_t());
		recycleQueue.enqueue_bulk(std::move_iterator<iter_t>(buffer.begin()),
				buffer.size());
	}

#pragma omp parallel
	{
		std::vector<kseq_t> readBuffer(bulkSize);
		string outBuffer;
		if (omp_get_thread_num() < filenames.size()) {
			//file reading init
			gzFile fp;
			fp = gzopen(filenames.at(omp_get_thread_num() - 1).c_str(), "r");
			kseq_t *seq = kseq_init(fp);

			//per thread token
			moodycamel::ProducerToken ptok(workQueue);

			//tokens for recycle queue
			moodycamel::ConsumerToken rctok(recycleQueue);

			unsigned dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
					std::move_iterator<iter_t>(readBuffer.begin()), bulkSize);
			while (dequeueSize == 0) {
				dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
						std::move_iterator<iter_t>(readBuffer.begin()),
						bulkSize);
			}

			unsigned size = 0;
			while (kseq_read(seq) >= 0) {
#pragma omp atomic update
				readCount++;
				cpy_kseq(&readBuffer[size++], seq);
				if (dequeueSize == size) {
					//try to insert, if cannot queue is full
					while (!workQueue.try_enqueue_bulk(ptok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							size)) {
						//try to work
						if (kseq_read(seq) >= 0) {
#pragma omp atomic update
							readCount++;
							consumer.consume(*seq);
#pragma omp atomic update
							processedReadsCount++;
						} else {
							goto fileEmpty;
						}
					}
					//reset buffer
					dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							bulkSize);
					while (dequeueSize == 0) {
						//try to work
						if (kseq_read(seq) >= 0) {
#pragma omp atomic update
							readCount++;
							consumer.consume(*seq);
#pragma omp atomic update
							processedReadsCount++;
						} else {
							goto fileEmpty;
						}
						dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
								std::move_iterator<iter_t>(readBuffer.begin()),
								bulkSize);
					}
					size = 0;
				}
			}
			fileEmpty:
			//finish off remaining work
			for (unsigned i = 0; i < size; ++i) {
				consumer.consume(readBuffer[i]);
			}
#pragma omp atomic update
			processedReadsCount += size;
			moodycamel::ProducerToken rptok(recycleQueue);
			assert(
					recycleQueue.enqueue_bulk(rptok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							size));
#pragma omp atomic update
			finished++;
			kseq_destroy(seq);
			gzclose(fp);
			//join in if others are still not finished
			if (finished < filenames.size()
					|| processedReadsCount < readCount) {
				moodycamel::ConsumerToken ctok(workQueue);
				while (finished < filenames.size()
						|| processedReadsCount < readCount) {
					size_t num = workQueue.try_dequeue_bulk(ctok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							bulkSize);
					if (num) {
						for (unsigned i = 0; i < num; ++i) {
							consumer.consume(readBuffer[i]);
						}
#pragma omp atomic update
						processedReadsCount += num;
						assert(
								recycleQueue.enqueue_bulk(rptok,
										std::move_iterator<iter_t>(
												readBuffer.begin()), num));
					}
				}
			}
		} else {
			moodycamel::ConsumerToken ctok(workQueue);
			moodycamel::ProducerToken rptok(recycleQueue);
			while (finished < filenames.size()
					|| processedReadsCount < readCount) {
				if (workQueue.size_approx() >= bulkSize) {
					size_t num = workQueue.try_dequeue_bulk(ctok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							bulkSize);
					if (num) {
						for (unsigned i = 0; i < num; ++i) {
							consumer.consume(readBuffer[i]);
						}
#pragma omp atomic update
						processedReadsCount += num;
						assert(
								recycleQueue.enqueue_bulk(rptok,
										std::move_iterator<iter_t>(
												readBuffer.begin()), num));
					}
				}
			}
		}
	}
}

#endif /* PRODCONKSEQRUNNER_HPP_ */

