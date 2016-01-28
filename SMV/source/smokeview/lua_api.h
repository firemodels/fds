#include "lua.h"

int lua_render(lua_State  *L);
void initLua();

int loadLuaScript();
void runLuaScript();
int loadSSFScript();
void runSSFScript();

void runScriptString(char *string);
int lua_get_sliceinfo(lua_State *L);
int lua_get_csvinfo(lua_State *L);
