% Calculate Gravitation load on the joints for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisOL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisOL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisOL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:58
% EndTime: 2020-06-27 17:29:59
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->9), mult. (26->14), div. (0->0), fcn. (16->4), ass. (0->6)
t12 = cos(qJ(1));
t11 = cos(qJ(3));
t10 = sin(qJ(1));
t9 = sin(qJ(3));
t8 = g(1) * t12 + g(2) * t10;
t1 = [(-MDP(3) + MDP(4)) * t8 + (MDP(6) * qJ(2) + MDP(2) + MDP(5)) * (-g(1) * t10 + g(2) * t12); t8 * MDP(6); (-g(1) * t9 + g(2) * t11) * MDP(8) + (-g(1) * t11 - g(2) * t9) * MDP(9); 0;];
taug = t1;
