/*! \file key_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include <functional>

class KeyLogger {
    public:
        using Code = int;
        using Callback = std::function<void(Code)>;

        KeyLogger();
        ~KeyLogger();

        bool install(Callback callback);
        bool terminate();
        const Callback & getCallback() const;

        static const char * codeToText(const Code & code);

    private:
        struct Data;
        std::unique_ptr<Data> data_;
};
