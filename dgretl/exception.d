/* The exception that every failure in the dgretl modules is reported with; its message has the details. */

module dgretl.exception;

class GretlException: Exception {
    this(string msg, string file = __FILE__, size_t line = __LINE__) {
        super(msg, file, line);
    }
}
