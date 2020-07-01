% Calculate Gravitation load on the joints for
% fourbar2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar2DE1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:54
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar2DE1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(2,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'fourbar2DE1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [3x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:54:47
% EndTime: 2020-06-26 17:54:47
% DurationCPUTime: 0.03s
% Computational Cost: add. (3->3), mult. (6->6), div. (0->0), fcn. (4->2), ass. (0->3)
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t1 = [(g(1) * t3 - g(2) * t4) * MDP(2) + (g(1) * t4 + g(2) * t3) * MDP(3);];
taug = t1;
