# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


class MetaboliteNotFound(Exception):
    """Thrown when a search for a given metabolite identifier does not yield any result."""
    pass


class PartNotFound(Exception):
    """Thrown by the ICE client when a requested genotype part is not found"""
    pass


class ModelStorageError(IOError):
    """Base exception thrown when a model can not be retrieved from the storage."""
    def __init__(self, message, *args, **kwargs):
        super().__init__(message, *args, **kwargs)
        self.message = message


class ModelNotFound(ModelStorageError):
    """Thrown when requesting a model which is not found."""
    pass


class Unauthorized(ModelStorageError):
    """Thrown when requesting a private model with invalid credentials."""
    pass


class Forbidden(ModelStorageError):
    """Thrown when requesting a private model for which the provided credentials are not authorized."""
    pass
