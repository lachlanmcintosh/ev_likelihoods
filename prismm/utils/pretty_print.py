def pretty_print(obj, max_length=1000):
    text = str(obj)
    if len(text) > max_length:
        text = text[:max_length] + '...'
    print(text)
