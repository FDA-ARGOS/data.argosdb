# code from: https://stackoverflow.com/a/60425679/2644759
import click
import questionary

from typing import Any, Optional, Sequence


class QuestionaryOption(click.Option):

    def __init__(
        self,
        param_decls: Optional[Sequence[str]] = None,
        **attrs: Any
    ) -> None:
        click.Option.__init__(self, param_decls, **attrs)
        if not isinstance(self.type, click.Choice):
            raise Exception('ChoiceOption type arg must be click.Choice')

    def prompt_for_value(self, ctx: click.Context) -> Any:
        prompt: str = self.prompt or 'Select one'

        if not isinstance(self.type, click.Choice):
            raise Exception('ChoiceOption type arg must be click.Choice')

        val = (
            questionary
            .select(prompt, choices=self.type.choices)
            .unsafe_ask()
        )
        return val
