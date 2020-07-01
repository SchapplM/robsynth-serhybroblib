% Calculate Gravitation load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnOL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnOL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:29
% EndTime: 2020-06-27 16:56:30
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->20), mult. (84->30), div. (0->0), fcn. (68->8), ass. (0->12)
t18 = qJ(2) + qJ(3);
t16 = sin(t18);
t17 = cos(t18);
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t26 = g(1) * t24 + g(2) * t21;
t27 = (g(3) * t17 - t26 * t16) * MDP(16) + (-g(3) * t16 - t26 * t17) * MDP(17);
t23 = cos(qJ(2));
t22 = cos(qJ(4));
t20 = sin(qJ(2));
t19 = sin(qJ(4));
t1 = [t26 * MDP(3) + (-t20 * MDP(10) - MDP(16) * t17 + MDP(17) * t16 + t22 * MDP(23) - t19 * MDP(24) + t23 * MDP(9) + MDP(2)) * (g(1) * t21 - g(2) * t24); (-g(3) * t23 + t26 * t20) * MDP(9) + (g(3) * t20 + t26 * t23) * MDP(10) + t27; t27; (-g(3) * t22 + t26 * t19) * MDP(23) + (g(3) * t19 + t26 * t22) * MDP(24); 0;];
taug = t1;
