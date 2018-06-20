#include <iostream>
#include <seqan/arg_parse.h>

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("modify_string");

    // set some documentation
    seqan::setShortDescription(parser, "String Modifier");
    seqan::setVersion(parser, "1.0");
    seqan::setDate(parser, "July 2012");

    /*
        Formatting Command Line Documentation

        The formatting of command line parameters might seem strange, at first: Font operators start with \f (which means that they start with "\\f" in in C++ string literals). The \\f is followed by the format specifier. The format specifier can be one of I, B, and P. I selects italic text (underlined on the shell), B selects bold and P resets the formatting to normal text. These font operators are legacies of man pages from Unix and offered a simple-to-implement solution to text formatting.

        For example, "Words \\fBwere\\fP made for \\fIbeing\\fP written!" would result in the formatted string “Words were made for being written!”.

        Note that formatting the command line relies on ANSI escape codes which is not supported by modern Windows versions. If you are using Windows, you will not see bold or underlined text.
    */
    seqan::addUsageLine(parser,
        "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    seqan::addDescription(parser,
        "This program allows simple character modifications to "
        "each i-th character.");

    seqan::addTextSection(parser, "Examples");

    seqan::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-U\\fP \\fIveryverylongword\\fP",
        "Print upper case version of \"veryverylongword\"");
    seqan::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP \\fIveryverylongword\\fP",
        "Print \"veryverylongword\" with every third character "
        "converted to upper case.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "period", 2);

    addOption(parser, seqan::ArgParseOption(
        "U", "uppercase", "Select to-uppercase as operation."));
    
    // this is a list
    addOption(parser, seqan::ArgParseOption(
        "A", "a", "show as list.", seqan::ArgParseArgument::INTEGER, "show as list", true));

    seqan::addOption(parser, seqan::ArgParseOption(
        "Z", "integer-value", "An integer option",
        seqan::ArgParseArgument::INTEGER, "INT"));

    // min/max values
    seqan::setMinValue(parser, "i", "10");
    seqan::setMaxValue(parser, "integer-value", "20");

    // required argument
    seqan::addOption(parser, seqan::ArgParseOption(
       "ir", "required-integer", "An required integer option",
       seqan::ArgParseArgument::INTEGER, "INT"));

   	setRequired(parser, "ir");

    // define list of valid values
    seqan::addOption(parser, seqan::ArgParseOption(
        "", "distance-model", "Distance model, either HAMMING or EDIT.",
    seqan::ArgParseArgument::STRING, "STR"));
    seqan::setValidValues(parser, "distance-model", "HAMMING EDIT");

    // define range of values
    seqan::addOption(parser, seqan::ArgParseOption(
        "r", "range", "The range to modify.",
        seqan::ArgParseArgument::INTEGER, "BEGIN END",
        false, 2)
    );

    // defining custom input/output file names
    seqan::addOption(parser, seqan::ArgParseOption(
        "I", "input-file", "Path to the input file",
        seqan::ArgParseArgument::INPUT_FILE, "IN"));
    seqan::addOption(parser, seqan::ArgParseOption(
        "O", "output-file", "Path to the output file",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    // allowed file ext
    seqan::setValidValues(parser, "input-file", "txt");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    unsigned period = 0;
    getOptionValue(period, parser, "period");
    bool toUppercase = isSet(parser, "uppercase");
    std::string text;
    int data;
    getArgumentValue(text, parser, 0);
    getOptionValue(data, parser, "a", 0);

    // get the range value
    unsigned rangeBegin = 0, rangeEnd = 0;
    seqan::getOptionValue(rangeBegin, parser, "range", 0);
    seqan::getOptionValue(rangeEnd, parser, "range", 1);

    std::cout << "period   \t" << period << '\n'
              << "uppercase\t" << toUppercase << '\n'
              << "text     \t" << text << '\n'
              << "data     \t" << data << '\n';

    return 0;
}