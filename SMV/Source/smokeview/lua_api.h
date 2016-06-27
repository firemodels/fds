#include "lua.h"

int lua_render(lua_State  *L);
void initLua();

int load_script(char *filename);
int loadLuaScript(char *filename);
int runLuaScript();
int loadSSFScript(char *filename);
int runSSFScript();

void runScriptString(char *string);
int lua_get_sliceinfo(lua_State *L);
int lua_get_csvinfo(lua_State *L);
