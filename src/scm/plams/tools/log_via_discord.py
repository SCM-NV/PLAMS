"""
resources: https://www.youtube.com/watch?v=CHbN_gB30Tw&list=PL-7Dfw57ZZVQ-GCNQS4Kyz637Fffhb0Hs&index=1&ab_channel=JamesS
"""

import json
from dataclasses import asdict, dataclass
from typing import Optional, Union

from scm.plams.core.functions import requires_optional_package


@dataclass
class SettingsDiscordLogger:
    """
    last update 2025
    channel_id:
        To get the channel_id go to advanced settings and toggle on the Dev mode, right click on a channel and the you can see the id

    token: https://www.youtube.com/watch?v=CHbN_gB30Tw&list=PL-7Dfw57ZZVQ-GCNQS4Kyz637Fffhb0Hs&index=1&ab_channel=JamesS
        search for `discord developer portal` on google
        log in, new application,
        switch on in page Bot the Presence intent, sever intent and message content intent
        reset token: will provide the token string to be copied
    """

    channel_id: int
    token: str

    def to_json(self, path="sett.json"):
        with open(path, "w") as f:
            json.dump(asdict(self), f)

    @classmethod
    def from_json(cls, path="sett.json"):
        with open(path, "r") as f:
            return cls(**json.load(f))

    def log_message(self, message: str):
        import discord

        class DiscordClient(discord.Client):

            def __init__(
                self,
                channel_id: int,
                startup_message: Optional[str] = None,
                file_paths: Optional[str] = None,
                close_on_start: bool = False,
                *,
                intents: discord.Intents,
                **options,
            ) -> None:
                self.file_paths = file_paths
                self.channel_id = channel_id
                self.startup_message = startup_message or f"🤖 Bot online!"
                self.close_on_start = close_on_start

                super().__init__(intents=intents, **options)

            async def on_ready(self):
                print(f"Logged on as {self.user}")

                channel = None
                if self.channel_id is not None:
                    channel = self.get_channel(self.channel_id)
                if channel is None:
                    print(f"❌ Could not find a channel with ID {self.channel_id}")
                else:
                    await channel.send(self.startup_message)
                    print(f"✅ Sent startup message in #{channel.name}, {self.startup_message}")

                if self.close_on_start:
                    await self.close()

        intents = discord.Intents.default()
        intents.message_content = True
        client = DiscordClient(
            channel_id=self.channel_id,
            intents=intents,
            close_on_start=True,
            startup_message=message,
        )
        client.run(self.token)
        return


@requires_optional_package("discord-py")
def log_message_via_discord(settings_path: Union[str, SettingsDiscordLogger], message: str):
    """_summary_

    :param settings_path: if is a str it should be a path to a json file with dict entries equivalent to attributes of SettingsDiscordLogger
    :type settings_path: Union[str, SettingsDiscordLogger]
    :param message: a string with a message to log
    :type message: str
    """
    if isinstance(settings_path, SettingsDiscordLogger):
        slog = settings_path
    else:
        slog = SettingsDiscordLogger.from_json(settings_path)
    slog.log_message(message)
