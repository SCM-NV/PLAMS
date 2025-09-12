"""
resources: https://www.youtube.com/watch?v=CHbN_gB30Tw&list=PL-7Dfw57ZZVQ-GCNQS4Kyz637Fffhb0Hs&index=1&ab_channel=JamesS
"""

import argparse
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import discord
import requests


def read_tail(p, n_lines=1):
    return os.popen(f"tail -n {n_lines} {p}").read()


class ExperimentLogger(discord.Client):
    options_keys = ("log-$int", "check", "log", "paths", "exit", "kill-all", "echo")

    def __init__(
        self,
        startup_message: Optional[str] = None,
        channel_id: Optional[int] = None,
        file_paths: Optional[str] = None,
        n_lines: int = 1,
        close_on_start: bool = False,
        message_tag: Optional[str] = None,
        *,
        intents: discord.Intents,
        **options,
    ) -> None:
        self.file_paths = file_paths
        self.n_lines = n_lines
        self.channel_id = channel_id
        self.message_tag = message_tag
        self.startup_message = (
            startup_message
            or f"🤖 Bot online! To activate use `{self.message_tag}` follow by options:`{'|'.join(self.options_keys)}`"
        )
        self.close_on_start = close_on_start

        super().__init__(intents=intents, **options)

    @property
    def all_path(self):
        if self.file_paths is None:
            return None
        return Path(self.file_paths).read_text().split("\n")

    def get_table_log(self, get_logs=False, n_lines=None):
        if self.all_path is None:
            return "Path are None"

        yield f"N Paths to log {len(self.all_path)}"
        for i, p in enumerate(self.all_path):
            ret = []
            ret += ["".join(["&"] * 40) + f" |Path {i}| " + "".join(["&"] * 40)]
            ret += [f"{p}"]
            ret += ["".join(["."] * 40) + "Tail" + "".join(["."] * 40)]
            if get_logs:
                ret += [f"{read_tail(p, n_lines or self.n_lines)}"]
            yield "```" + "\n".join(ret) + "```"
        # ret = {"Paths": self.all_path}
        # maxcolwidths = [20]
        # if get_logs:
        #     ret["Logs"] = [read_tail(p, n_lines or self.n_lines) for p in self.all_path]
        #     maxcolwidths.append(70)
        # table = tabulate(ret, headers="keys",  maxcolwidths=maxcolwidths)
        # return f"```\n{table}\n```"

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

    async def on_message(self, message):
        # don't respond to ourselves
        if message.author == self.user:
            return

        if message.content == "kill-all":
            await self.close()

        if message.content == "echo":
            await message.channel.send(self.message_tag or "No Tag")

        if self.message_tag is not None:
            if self.message_tag not in message.content:
                print(f"{self.message_tag=} not in {message.content}")
                return
            message.content = message.content.replace(self.message_tag, "").strip()

        if "log-" in message.content:
            n_lines = int(message.content.split("-")[-1])
            for m in self.get_table_log(get_logs=True, n_lines=n_lines):
                await message.channel.send(m)

        if message.content in ["check", "log"]:
            for m in self.get_table_log(get_logs=True):
                await message.channel.send(m)

        if message.content == "paths":
            for m in self.get_table_log():
                await message.channel.send(m)

        if message.content == "exit":
            await self.close()


@dataclass
class SettingsLogger:
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

    def send_message(
        self,
        auth: str,
        message: str = "Hi",
    ):
        # source https://www.youtube.com/watch?v=g0jsaTIWz5I&ab_channel=CodeBear
        # auth is not self.token

        url = f"https://discord.com/api/v9/channels/{self.channel_id}/messages"
        payload = {"content": message}
        headers = {"Authorization": auth}
        res = requests.post(url, payload, headers=headers)
        print("sent something")


def log_message(message: str, settings_path: Optional[str] = None):
    if settings_path is not None:
        slog = SettingsLogger.from_json(settings_path)
    else:
        slog = SettingsLogger.from_json()
    intents = discord.Intents.default()
    intents.message_content = True
    client = ExperimentLogger(
        channel_id=slog.channel_id,
        intents=intents,
        close_on_start=True,
        startup_message=message,
    )
    client.run(slog.token)
    return


def app_main(
    settings_path: Optional[str] = None, file_paths="ExperimentLogger_paths.txt", message_tag: Optional[str] = None
):
    """Launch an ap that check tail of every file listed in ExperimentLogger_paths.txt"""
    if settings_path is not None:
        slog = SettingsLogger.from_json(settings_path)
    else:
        slog = SettingsLogger.from_json()

    intents = discord.Intents.default()
    intents.message_content = True
    client = ExperimentLogger(
        channel_id=slog.channel_id,
        file_paths=file_paths,
        intents=intents,
        message_tag=message_tag,
    )
    client.run(slog.token)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--settings-path", default=None, help="Path to a json file with the settings")
    parser.add_argument(
        "-e",
        "--experiments-path",
        default="ExperimentLogger_paths.txt",
        help="Path to txt file with each line a path to a file to be monitored",
    )
    parser.add_argument(
        "-t",
        "--message-tag",
        default=None,
        help="Path to txt file with each line a path to a file to be monitored",
    )
    args = parser.parse_args()
    print(args.settings_path)
    app_main(settings_path=args.settings_path, file_paths=args.experiments_path, message_tag=args.message_tag)


if __name__ == "__main__":
    main()
