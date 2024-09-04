def escape_text(text):
    return text.translate(str.maketrans(
        {"\"":'""',
        # "'": "\\'",
        "\\": "\\\\",
        ";": "\\;"
        }
        ))