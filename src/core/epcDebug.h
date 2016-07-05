#ifndef EPC_DEBUG_H
#define EPC_DEBUG_H

#include <cassert>

#ifdef EPC_DEBUG
    #define EPC_ASSERT( COND )        assert( COND );

	#define EPC_MSG_ASSERT(msg, cond) do \
	{ if (!(cond)) { std::ostringstream str; str << msg; std::cerr << str.str(); std::abort(); } \
	} while(0)
#else
    #define EPC_ASSERT( COND )

    #define EPC_MSG_ASSERT(msg, cond)
#endif  // EPC_DEBUG

#endif // EPC_DEBUG_H
