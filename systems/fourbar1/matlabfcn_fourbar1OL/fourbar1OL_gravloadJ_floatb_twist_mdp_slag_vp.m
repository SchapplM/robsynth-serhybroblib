% Calculate Gravitation load on the joints for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1OL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1OL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1OL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:19
% EndTime: 2020-06-26 17:43:19
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->11), mult. (24->18), div. (0->0), fcn. (16->6), ass. (0->9)
t14 = qJ(1) + qJ(2);
t12 = sin(t14);
t13 = cos(t14);
t19 = (-g(1) * t12 + g(2) * t13) * MDP(5) + (-g(1) * t13 - g(2) * t12) * MDP(6);
t18 = cos(qJ(1));
t17 = cos(qJ(3));
t16 = sin(qJ(1));
t15 = sin(qJ(3));
t1 = [(g(1) * t16 - g(2) * t18) * MDP(2) + (g(1) * t18 + g(2) * t16) * MDP(3) + t19; t19; (g(1) * t15 - g(2) * t17) * MDP(8) + (g(1) * t17 + g(2) * t15) * MDP(9); 0;];
taug = t1;
