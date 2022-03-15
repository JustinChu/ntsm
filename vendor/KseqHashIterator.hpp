/*
 * KseqHash.hpp
 *
 *  Created on: Sep 28, 2021
 *      Author: cjustin
 */

#ifndef KSEQHASHITERATOR_HPP_
#define KSEQHASHITERATOR_HPP_

#include <stdint.h>
#include <string>
#include <limits>

class KseqHashIterator {
public:

	/**
	 * Default constructor. Creates an iterator pointing to
	 * the end of the iterator range.
	 */
	KseqHashIterator() :
			m_seq(NULL), m_len(0), m_k(0), m_mask(0), m_shift(0), m_pos(
					std::numeric_limits<std::uint64_t>::max()) {
		init();
	}

	KseqHashIterator(const char *seq, uint64_t strLen, int k) :
			m_seq(seq), m_len(strLen), m_k(k), m_mask((1ULL << k * 2) - 1), m_shift(
					(k - 1) * 2), m_pos(0) {
		init();
		next();
	}

	/** get pointer to hash value for current k-mer */
	uint64_t operator*() const {
		return m_hashVal;
	}

	/** test equality with another iterator */
	bool operator==(const KseqHashIterator &it) const {
		return m_pos == it.m_pos;
	}

	/** test inequality with another iterator */
	bool operator!=(const KseqHashIterator &it) const {
		return !(*this == it);
	}

	/** pre-increment operator */
	KseqHashIterator& operator++() {
		next();
		return *this;
	}

	static const KseqHashIterator end() {
		return KseqHashIterator();
	}

	uint64_t getPos(){
		return m_pos;
	}

	~KseqHashIterator() {
	}

private:
	const char *m_seq;
	const uint64_t m_len;
	const unsigned m_k;
	const uint64_t m_mask;
	const uint64_t m_shift;
	uint64_t m_pos;
	unsigned m_subStrLen;
	uint64_t m_ntFW;
	uint64_t m_ntRV;
	uint64_t m_hashVal;

	/** Initialize internal state of iterator */
	void init() {
		m_subStrLen = 0;
		m_ntFW = 0;
		m_ntRV = 0;
		m_hashVal = 0;
	}

	void next() {
		step();
		while (m_subStrLen < m_k
				&& m_pos != std::numeric_limits<uint64_t>::max()) {
			step();
		}
	}

	void step() {
		if (m_pos < m_len) {
			int c = s_seq_nt4_table[(uint8_t) m_seq[m_pos++]];
			if (c < 4) { // not an "N" base
				m_ntFW = (m_ntFW << 2 | c) & m_mask;       // forward strand
				m_ntRV = m_ntRV >> 2 | (uint64_t) (3 - c) << m_shift; // reverse strand
				if (++m_subStrLen >= m_k) { // we find a k-mer
					m_hashVal = hash64(m_ntFW < m_ntRV ? m_ntFW : m_ntRV,
							m_mask);
					return;
				}
			} else
				init(); // if there is an "N", restart
		}
		else{
			m_pos = std::numeric_limits<uint64_t>::max();
		}
	}

	const unsigned char s_seq_nt4_table[256] = { // translate ACGT to 0123
			0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4,
					2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

	inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
			{
		key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
		key = key ^ key >> 24;
		key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
		key = key ^ key >> 14;
		key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
		key = key ^ key >> 28;
		key = (key + (key << 31)) & mask;
		return key;
	}
};

#endif /* KSEQHASHITERATOR_HPP_ */
