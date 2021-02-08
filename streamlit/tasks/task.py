

class Task():
    def __init__(self, render_kwargs, compute_kwargs, next_task):
        self.kwargs = render_kwargs

    def run(self):
        self.render()
        self.compute()

    def render(self):
        pass

    def compute(self):
        pass
