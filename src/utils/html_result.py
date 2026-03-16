import sys
from utils.logger_debug import setup_log


class HTMLResult(object):
    def __init__(self, output="", as_dev=False):
        self.logger = setup_log(__name__, as_dev)
        self.output_name = output
        self.header = ""
        self.footer = f'</body> </html>'
        self.content = ""
        self.css = ""
        self.js = ""

    def make_page(self):
        return "\n".join([self.header, self.css, self.js, self.content, self.footer])

    def make_header(self):
        self.header = f'<html> <head><title>Results {self.output_name}</title> {self.css} {self.css}  </head> <body>'

    def set_content(self, html=""):
        self.content = html


