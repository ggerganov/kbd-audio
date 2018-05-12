/*! \file key_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include <functional>

class KeyLogger {
    public:
        using Code = int;
        using Callback = std::function<void(Code)>;

        KeyLogger();

        bool install(Callback callback);
        const Callback & getCallback() const;

        static const char * codeToText(const Code & code);

    private:
        Callback callback_;
};
